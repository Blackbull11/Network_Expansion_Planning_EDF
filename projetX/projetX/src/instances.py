"""
Lecture des instances, quel que soit leur format :
  - data/instances/*.json (generate_instances.py) : p_max_depart et
    l_max_depart sur chaque source, demande_hivers présent ;
  - verified_instances/*.json : limites globales dans le bloc
    'electrotechnique' (puissance_max_ligne, longueur_max_depart), pas de
    demande_hivers (un fichier par saison).

load_instance renvoie un GraphData utilisable tel quel par voronoi.py,
steiner.py et verify.py : si le fichier n'a pas le format de
data/instances, une copie normalisée est écrite dans NORMALIZED_DIR.
"""
import json
from pathlib import Path

from nx_graph import GraphData

ROOT = Path(__file__).resolve().parent.parent
VERIFIED_DIR = ROOT / "verified_instances"
NORMALIZED_DIR = ROOT / "output" / "normalized_instances"
SIZE_CLASSES = ('Petit', 'Moyen', 'Grand', 'TrèsGrand')   # first word of the names in verified_instances/


def normalize(data: dict) -> bool:
    """
    Met `data` (JSON d'une instance) au format de data/instances, en place.
    Renvoie True si quelque chose a changé.
    """
    elec = data.get('electrotechnique', {})
    defaults = {'p_max_depart': elec.get('puissance_max_ligne'),
                'l_max_depart': elec.get('longueur_max_depart')}
    changed = False
    for n in data['nodes']:
        if 'demande_hivers' not in n:
            n['demande_hivers'] = n['demande']
            changed = True
        if n['type'] != 'Source':
            continue
        for key, value in defaults.items():
            if n.get(key) is None and value is not None:
                n[key] = value
                changed = True
    return changed


def load_instance(path) -> GraphData:
    """GraphData de l'instance `path`, normalisée si besoin (voir normalize)."""
    path = Path(path)
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    if normalize(data):
        NORMALIZED_DIR.mkdir(parents=True, exist_ok=True)
        path = NORMALIZED_DIR / path.name
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data, f)
    return GraphData(folderpath=str(path.parent), filename=path.name)


def instance_size(path) -> tuple:
    """(nb de nœuds, nb d'arêtes, nb de sources) lus dans le JSON, pour trier les instances."""
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    n_src = sum(n['type'] == 'Source' for n in data['nodes'])
    return len(data['nodes']), len(data['edges']), n_src


def is_instance(path) -> bool:
    """True si `path` est bien un JSON d'instance (certains .json de verified_instances/ sont des PDF ou des CSV)."""
    try:
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
    except (UnicodeDecodeError, json.JSONDecodeError):
        return False
    return isinstance(data, dict) and 'nodes' in data and 'edges' in data


def verified_files(prefix: str = "") -> list:
    """Instances JSON de verified_instances/ dont le nom commence par `prefix`, triées par taille réelle."""
    files = [f for f in VERIFIED_DIR.glob(f"{prefix}*.json") if f.stem.split('_')[0] in
             SIZE_CLASSES and is_instance(f)]
    return sorted(files, key=lambda f: (instance_size(f), f.name))
