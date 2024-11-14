import json
import logging
from http.client import HTTPResponse
from typing import Optional, List

from fastapi import APIRouter, HTTPException, UploadFile, File, Request
from fastapi.openapi.models import Response
from rdkit import Chem

from app.config import settings
from app.models.molecules import Molecule
from app.serializers.serializers import MoleculeCreateSerializer, MoleculeReadSerializer, validate_smiles, \
    MoleculeUpdateSerializer, MoleculeListCreateSerializer, MoleculeListSerializer
from app.storage.cache import get_cached_result, set_cache
from app.storage.postgres.molecules import MoleculeStorage
import utils.logger

logger = logging.getLogger(__name__)

mol_router = APIRouter()

@mol_router.get("/")
def get_test_message():
    logger.info('get_test_message...')
    return {"instance_name": settings.INSTANCE_NAME}


@mol_router.post("/molecules", response_model=MoleculeReadSerializer)
def create_molecule(molecule_data: MoleculeCreateSerializer):
    logger.info(f'create_molecule...\n{molecule_data}')
    storage = MoleculeStorage()
    try:
        molecule = Molecule(**molecule_data.model_dump())
        molecule = storage.save(molecule)
        return molecule
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@mol_router.get("/molecules/{molecule_id}", response_model=MoleculeReadSerializer)
def get_molecule_by_id(molecule_id: int):
    logger.info(f'get_molecule_by_id...\nid: {molecule_id}')
    storage = MoleculeStorage()
    molecule = storage.get_by_id(molecule_id)
    if molecule:
        return molecule
    else:
        raise HTTPException(status_code=404, detail="Molecule not found")


# @mol_router.get("/molecules", response_model=MoleculeReadSerializer)
# def get_molecule_by_smiles(smiles: str):
#     try:
#         validate_smiles(smiles)
#     except ValueError as e:
#         raise HTTPException(status_code=400, detail=str(e))
#
#     storage = MoleculeStorage()
#     molecule = storage.get_by_smiles(smiles)
#     if molecule:
#         return Response(status=200, content=molecule)
#     else:
#         raise HTTPException(status_code=404, detail="Molecule not found")


@mol_router.put("/molecules/{molecule_id}", response_model=MoleculeReadSerializer)
def update_molecule(molecule_id: int, molecule_data: MoleculeUpdateSerializer):
    logger.info(f'update_molecule...\nmolecule_id: {molecule_id}, molecule_data: {molecule_data}')
    storage = MoleculeStorage()
    try:
        molecule = Molecule(**molecule_data.model_dump())
        molecule = storage.partial_update(molecule_id, molecule)
        return molecule
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@mol_router.delete("/molecules/{molecule_id}", response_model=dict)
def delete_molecule(molecule_id: int):
    logger.info(f'delete_molecule...\nmolecule_id: {molecule_id}')
    storage = MoleculeStorage()
    try:
        storage.delete(molecule_id)
        return {"detail": f"Molecule with ID {molecule_id} deleted"}
    except ValueError as e:
        logger.warning(f'delete_molecule {str(e)}')
        raise HTTPException(status_code=404, detail=str(e))


@mol_router.get("/molecules", response_model=List[MoleculeReadSerializer])
def search_molecules(request: Request, smiles: Optional[str] = None, limit: Optional[int] = 1000, offset: Optional[int] = 0):
    logger.info(f'search_molecules...\nsmiles: {smiles}, limit: {limit}, offset: {offset}')

    cache_key = f'search:{request.url}'
    cached_result = get_cached_result(cache_key)
    if cached_result:
        logger.info(f'cached_result: {cached_result}')
        return cached_result

    storage = MoleculeStorage()

    if smiles:
        try:
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
            set_cache(cache_key, [molecule.model_dump() for molecule in molecule_list_schema])
            set_cache(cache_key, [result.__dict__ for result in results])
            return results
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))

    try:
        results = storage.get_all()
        molecule_list_schema = [MoleculeReadSerializer.model_validate(result.__dict__) for result in results]
        set_cache(cache_key, [molecule.model_dump() for molecule in molecule_list_schema])
        return results
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))



@mol_router.post("/molecules/upload-json", response_model=List[MoleculeReadSerializer])
async def upload_molecules_from_json(file: UploadFile = File(...)):
    logger.info(f'search_molecules...\nfile: {file.file.name}')
    storage = MoleculeStorage()

    if file.content_type != 'application/json':
        raise HTTPException(status_code=400, detail="Invalid file type. Only JSON files are supported.")

    try:
        contents = await file.read()
        json_data = json.loads(contents.decode('utf-8'))

        try:
            logger.info(f'json_data: {json_data}')
            molecule_data = MoleculeListCreateSerializer(**json_data)
        except Exception as e:
            logger.warning(f'upload_molecules_from_json Invalid JSON structure: {str(e)}')
            raise HTTPException(status_code=400, detail=f"Invalid JSON structure: {str(e)}")

        created_molecules = []
        for molecule_info in molecule_data.molecules:
            logger.info(f'molecule_info: {molecule_info}')
            try:
                mol = Molecule(**molecule_info.model_dump())
                logger.info(f'mol : {mol}')
                if mol is None:
                    raise ValueError(f"Invalid molecule structure: [{molecule_info.name}, {molecule_info.smiles}]")
                stored_molecule = storage.save(mol)
                created_molecules.append(stored_molecule)
            except Exception as e:
                logger.warning(f'upload_molecules_from_json Invalid JSON structure: {str(e)}')
                raise HTTPException(status_code=400, detail=str(e))

        return created_molecules
    except Exception as e:
        logger.error(f'upload_molecules_from_json Error: {str(e)}')
        raise HTTPException(status_code=500, detail=f"An error occurred while processing the file: {str(e)}")
