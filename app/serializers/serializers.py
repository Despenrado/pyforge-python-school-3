from pydantic import BaseModel, validator, ValidationError
from typing import Optional, List, TypeVar, Generic

from rdkit import Chem


def validate_smiles(value):
    if not Chem.MolFromSmiles(value):
        raise ValidationError('Invalid SMILES string')
    return value


class MoleculeCreateSerializer(BaseModel):
    name: str
    smiles: str

    @validator('smiles')
    def validate_smiles(cls, value):
        return validate_smiles(value)

    class Config:
        orm_mode = True


class MoleculeReadSerializer(BaseModel):
    id: int
    name: str
    smiles: str

    class Config:
        orm_mode = True


class MoleculeListSerializer(BaseModel):
    molecules: List[MoleculeReadSerializer]

    class Config:
        orm_mode = True


class MoleculeListCreateSerializer(BaseModel):
    molecules: List[MoleculeCreateSerializer]

    class Config:
        orm_mode = True


class MoleculeUpdateSerializer(BaseModel):
    name: Optional[str] = None
    smiles: Optional[str] = None

    @validator('smiles')
    def validate_smiles(cls, value):
        return validate_smiles(value)

    class Config:
        orm_mode = True

T = TypeVar('T')
class CeleryResult(BaseModel, Generic[T]):
    task_id: str
    status: str
    result: Optional[T] = None

    class Config:
        orm_mode = True