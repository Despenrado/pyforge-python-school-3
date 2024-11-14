from app.storage.postgres.db_connect import Base
from sqlalchemy import Column, Integer, String


class Molecule(Base):
    __tablename__ = 'molecules'

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, index=True)
    smiles = Column(String, index=True)
    # TODO: fingerprint

