
import logging
from app.celery_app import celery
from rdkit import Chem
from app.serializers.serializers import validate_smiles, MoleculeReadSerializer
from app.storage.cache import get_cached_result, set_cache
from app.storage.postgres.molecules import MoleculeStorage
import utils.logger

logger = logging.getLogger(__name__)

@celery.task
def perform_substructure_search(request_url, smiles, offset=0, limit=1000):
    cache_key = f'search:{request_url}'
    cached_result = get_cached_result(cache_key)
    if cached_result:
        logger.info(f'cached_result: {cached_result}')
        return cached_result

    storage = MoleculeStorage()

    if smiles:
        validate_smiles(smiles)
        substructure_smiles = Chem.MolFromSmiles(smiles)

        results = []
        for molecules_batch in storage.get_all(offset=offset, limit=limit):
            for molecule in molecules_batch:
                mol = Chem.MolFromSmiles(molecule.smiles)
                if mol and mol.HasSubstructMatch(substructure_smiles):
                    results.append(molecule)
        logger.info(f'search_molecules: {results}')
        molecule_list_schema = [MoleculeReadSerializer.model_validate(result.__dict__) for result in results]
        res = [molecule.model_dump() for molecule in molecule_list_schema]
        set_cache(cache_key, res)
        return res

    results = storage.get_all()
    molecule_list_schema = [MoleculeReadSerializer.model_validate(result.__dict__) for result in results]
    res = [molecule.model_dump() for molecule in molecule_list_schema]
    set_cache(cache_key, res)
    return res
