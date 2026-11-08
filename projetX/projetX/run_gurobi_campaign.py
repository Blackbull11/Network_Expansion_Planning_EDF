"""
Campagne Gurobi (plne_gurobi : modèle renforcé, arrêt à la borne de subset_dp)
en parallèle : 5 Petit_N75 et 5 Petit_N90 tirés au hasard, puis les Moyen des
plus petits aux plus grands, jusqu'à ce que TIMEOUT_TARGET runs Moyen atteignent
la limite de temps ou de mémoire sans optimum prouvé. Les runs en cours à ce
moment-là vont à leur terme.

Suivi : output/gurobi_campaign/progress.log. LAUNCH et FINISH à chaque run,
STATUS à chaque fin de run et toutes les 5 minutes, avec le temps restant estimé.
Autres résultats : runs/<instance>.json, runs/<instance>.gurobi.log, summary.csv.

Le script reprend là où il s'est arrêté (un run déjà dans runs/ n'est pas relancé).

Usage :
  python run_gurobi_campaign.py
  python run_gurobi_campaign.py --workers 3 --threads 5 --time-limit 1200 --seed 2026
"""
import argparse
import contextlib
import csv
import io
import json
import logging
import random
import time
import traceback
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
INSTANCE_DIR = HERE / "data" / "verified_instances"
OUT_DIR = HERE / "output" / "gurobi_campaign"
RUNS_DIR = OUT_DIR / "runs"
LOG = OUT_DIR / "progress.log"

N_PER_SIZE = 5
TIMEOUT_TARGET = 5
SOLVED = ('OPTIMAL', 'OPTIMAL_BY_DP')
TIMEOUT_STATUS = ('TIME_LIMIT', 'MEM_LIMIT')   # une solution non prouvée, faute de temps ou de mémoire
PARAMS = {'MIPGap': 1e-6, 'NodefileStart': 0.5, 'SoftMemLimit': 3.0}   # Go
SUMMARY_COLS = ('instance', 'family', 'n_nodes', 'n_sources', 'n_loads', 'dp_bound', 'status', 'obj', 'bound', 'gap',
                'feasible', 'time_dp', 'time_build', 'time_solve', 'time_total')


def log(message: str) -> None:
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {message}\n"
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line)
    print(line, end="", flush=True)


def hms(seconds: float) -> str:
    seconds = int(round(max(seconds, 0)))
    return f"{seconds // 3600:d}h{seconds % 3600 // 60:02d}m{seconds % 60:02d}s"


def load_count(path: Path) -> int:
    with open(path, encoding="utf-8") as fh:
        return sum(n['type'] in ('Conso', 'Prod') for n in json.load(fh)['nodes'])


def instance_list(seed: int) -> list:
    """[(stem, famille)] : les Petit tirés au hasard, puis tous les Moyen triés par nombre de charges."""
    rng = random.Random(seed)
    jobs = []
    for size in ("Petit_N75", "Petit_N90"):
        names = sorted(p.stem for p in INSTANCE_DIR.glob(f"{size}_*.json"))
        jobs += [(n, "Petit") for n in sorted(rng.sample(names, N_PER_SIZE))]
    moyen = sorted(INSTANCE_DIR.glob("Moyen_*.json"), key=lambda p: (load_count(p), p.name))
    return jobs + [(p.stem, "Moyen") for p in moyen]


def run_job(stem: str, family: str, a: dict) -> dict:
    """Un run complet : borne subset_dp si l'instance est assez petite, puis Gurobi et certificat."""
    t0 = time.time()
    logging.disable(logging.INFO)
    import plne_gurobi as P
    import subset_dp
    from src.instances import instance_size, load_instance
    from verify import check_solution
    path = INSTANCE_DIR / f"{stem}.json"
    with contextlib.redirect_stdout(io.StringIO()):
        graph = load_instance(path)
    n_loads = sum(lab in ('Conso', 'Prod') for _, lab in graph.G.nodes('label'))
    n_nodes, _, n_sources = instance_size(path)
    res = dict(instance=stem, family=family, n_nodes=n_nodes, n_sources=n_sources, n_loads=n_loads, dp_bound=None,
               time_dp=None, status=None, obj=None, bound=None, gap=None, feasible=None, time_build=None, time_solve=None)
    try:
        if n_loads <= P.DP_MAX_LOADS:
            dp = subset_dp.lower_bound(graph)
            res.update(dp_bound=dp['bound'], time_dp=dp['time'])
            log(f"DP      {stem}: lower bound {dp['bound']:.4f} km in {hms(dp['time'])}")
        if res['dp_bound'] == float('inf'):
            res.update(dp_bound=None, status='INFEASIBLE_BY_DP')
        else:
            solved, trees = P.solve_graph(graph, a['time_limit'], threads=a['threads'], warm=False,
                                          log_file=RUNS_DIR / f"{stem}.gurobi.log", tighten=True,
                                          dp_bound=res['dp_bound'], params=PARAMS)
            res.update(solved)
            if trees is not None:
                cert = check_solution(graph, trees)
                res.update(feasible=cert['feasible'], J_verify=cert['J'], errors=[str(e) for e in cert['errors']],
                           trees={str(r): es for r, es in trees.items()})
                P.certify_with_bound(res, cert)
    except Exception as exc:
        res.update(status=f"error: {exc!r}", traceback=traceback.format_exc())
    res['time_total'] = time.time() - t0
    with open(RUNS_DIR / f"{stem}.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, default=str)
    return res


class Campaign:

    def __init__(self, a):
        self.a = a
        self.total_petit = 2 * N_PER_SIZE
        self.queue = []            # [(stem, famille)] restant à lancer
        self.running = {}          # {future: (stem, famille, début)}
        self.rows = []
        self.durations = {'Petit': [], 'Moyen': []}
        self.moyen_done = 0
        self.timeouts = 0

    def mean(self, family):
        d = self.durations[family]
        return sum(d) / len(d) if d else (0.5 * self.a.time_limit if family == 'Petit' else float(self.a.time_limit))

    def estimate(self) -> tuple:
        """(temps restant estimé, pire cas) en secondes."""
        now = time.time()
        limit = self.a.time_limit
        q_petit = [j for j in self.queue if j[1] == 'Petit']
        q_moyen = [j for j in self.queue if j[1] == 'Moyen']
        work = sum(max(60.0, self.mean(f) - (now - t)) for _, f, t in self.running.values())
        worst = sum(max(60.0, limit - (now - t)) for _, _, t in self.running.values())
        work += len(q_petit) * self.mean('Petit')
        worst += (len(q_petit) + len(q_moyen)) * limit
        if self.timeouts < TIMEOUT_TARGET:
            rate = (self.timeouts + 1) / (self.moyen_done + 2)   # part de runs Moyen sans preuve observée
            needed = (TIMEOUT_TARGET - self.timeouts) / rate
            needed -= sum(1 for _, f, _ in self.running.values() if f == 'Moyen')
            work += min(max(needed, 0), len(q_moyen)) * self.mean('Moyen')
        else:
            worst -= len(q_moyen) * limit
        return work / self.a.workers, worst / self.a.workers

    def status_line(self) -> None:
        done_petit = sum(r['family'] == 'Petit' for r in self.rows)
        running = ", ".join(f"{s} ({hms(time.time() - t)})" for s, _, t in self.running.values()) or "none"
        expected, worst = self.estimate()
        log(f"STATUS  Petit {done_petit}/{self.total_petit} done | Moyen {self.moyen_done} done, "
            f"{self.timeouts}/{TIMEOUT_TARGET} unproven timeouts | queued {len(self.queue)}")
        log(f"        running: {running}")
        log(f"        estimated remaining: {hms(expected)} (worst case {hms(worst)})")

    def write_summary(self) -> None:
        with open(OUT_DIR / "summary.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, extrasaction='ignore')
            w.writeheader()
            w.writerows(self.rows)

    def record(self, res: dict) -> None:
        self.rows.append(res)
        if res['family'] == 'Moyen':
            self.moyen_done += 1
            self.timeouts += res['status'] in TIMEOUT_STATUS
        self.durations[res['family']].append(res.get('time_total') or 0.0)
        self.write_summary()


def main() -> None:
    p = argparse.ArgumentParser(description="Campagne Gurobi avec arrêt à la borne de subset_dp (voir le docstring)")
    p.add_argument("--workers", type=int, default=3, help="instances résolues en même temps")
    p.add_argument("--threads", type=int, default=5, help="threads Gurobi par instance")
    p.add_argument("--time-limit", type=float, default=1200, help="limite de temps Gurobi par instance (s)")
    p.add_argument("--seed", type=int, default=2026, help="graine du tirage des instances Petit")
    a = p.parse_args()
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    camp = Campaign(a)
    jobs = instance_list(a.seed)
    for stem, family in jobs:
        f = RUNS_DIR / f"{stem}.json"
        if f.exists():
            with open(f, encoding="utf-8") as fh:
                res = json.load(fh)
            if not str(res.get('status')).startswith('error'):
                camp.record(res)
                continue
        camp.queue.append((stem, family))
    log("=" * 88)
    log(f"CAMPAIGN START: {a.workers} runs at a time, {a.threads} threads each, {a.time_limit:.0f} s Gurobi limit per run.")
    log(f"Petit runs: {', '.join(s for s, f in jobs if f == 'Petit')}")
    log(f"Moyen runs: smallest first, until {TIMEOUT_TARGET} reach the time or memory limit without a proven optimum.")
    log("Moyen has 22 to 62 loads: the DP bound is only computed up to 20 loads, so no Moyen run has a DP bound.")
    log(f"Already done and skipped: {len(camp.rows)}.")
    args = dict(time_limit=a.time_limit, threads=a.threads)
    t_start = time.time()
    with ProcessPoolExecutor(max_workers=a.workers, max_tasks_per_child=1) as ex:
        while camp.queue or camp.running:
            while camp.queue and len(camp.running) < a.workers and camp.timeouts < TIMEOUT_TARGET:
                stem, family = camp.queue.pop(0)
                fut = ex.submit(run_job, stem, family, args)
                camp.running[fut] = (stem, family, time.time())
                log(f"LAUNCH  {stem} ({family}, {load_count(INSTANCE_DIR / f'{stem}.json')} loads)")
            if not camp.running:
                break
            done, _ = wait(list(camp.running), timeout=300, return_when=FIRST_COMPLETED)
            for fut in done:
                stem, family, t = camp.running.pop(fut)
                try:
                    res = fut.result()
                except Exception as exc:
                    res = dict(instance=stem, family=family, status=f"error: {exc!r}", time_total=time.time() - t)
                camp.record(res)
                lb = f"{res['dp_bound']:.4f}" if res.get('dp_bound') is not None else "none"
                obj = f"{res['obj']:.4f}" if res.get('obj') is not None else "none"
                solver_lb = f"{res['bound']:.4f}" if res.get('bound') is not None else "none"
                log(f"FINISH  {stem}: {res['status']} | J = {obj} km | DP bound {lb} | solver bound {solver_lb} | "
                    f"total {hms(res.get('time_total') or 0)}")
            camp.status_line()
    n_solved = sum(r['status'] in SOLVED for r in camp.rows)
    log(f"CAMPAIGN END: {len(camp.rows)} runs, {n_solved} proven optimal, {camp.timeouts} unproven Moyen timeouts, "
        f"in {hms(time.time() - t_start)}.")


if __name__ == "__main__":
    main()
