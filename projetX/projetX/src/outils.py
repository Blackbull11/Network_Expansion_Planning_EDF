"""
Fonctions utiles (et utilisées dans différents programmes).
"""

import json
import logging
import numpy as np
import pandas as pd
import geopandas

"""
Fonctions pour les SOLVEURS.
"""


def get_solution(var_decision: dict) -> dict:
    """
    Converti les variables linopy en leur solution - applicable une fois que le modèle a été résolu, s'il admet une solution.

    Args:
        - var_decision : dictionnaire contenant les variables du modèle résolu.

    Returns:
        - solution : un dictionnaire contenant les mêmes colonnes, où les variables linopy ont été remplacées par leur valeur
        pour la solution du problème.
    """
    solution = {}
    for vd_key, vd_array in var_decision.items():
        solution[vd_key] = vd_array.solution.to_numpy().tolist()
    return solution


def write_output(
    parameters: dict, feasible: bool, solution: dict, filename="simulation_output.json"
):
    """
    Ecrit la solution dans un fichier json.
    """
    output = {"parametres": parameters, "faisable": feasible, "solution": solution}
    with open(filename, "w") as f:
        json.dump(output, f)


"""
Fonctions pour les LOGS et la GESTION D'ERREUR.
"""


def init_stream_logger(level=logging.INFO) -> logging.Logger:
    """
    Initialise le logger pour avoir les logs dans le terminal.
    """
    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        level=level,
    )
    logger = logging.getLogger(__name__)
    #logger.addHandler(logging.StreamHandler())
    return logger


def add_file_logger(
    log_filename: str, level=logging.DEBUG, filemode="w"
) -> logging.Logger:
    """
    Ajoute un handler pour avoir les logs dans un fichier de logs.
    """
    logger = logging.getLogger(__name__)
    fh = logging.FileHandler(log_filename, mode=filemode)
    fh.setLevel(level)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger, fh


# def print_or_log_INFO(msg: str, sync_log: bool, logger: logging.Logger):
#     if sync_log:
#         logger.info(msg)
#     else:
#         print(msg)


def array_float_to_int_with_tol(my_array, tolerance=1e-6):
    arr_rounded = np.where(
        np.abs(my_array - np.round(my_array)) < tolerance, np.round(my_array), my_array
    )
    arr_int = arr_rounded.astype(int)
    return arr_int


if __name__ == "__main__":
    log = init_stream_logger()

    log.info("test info")
