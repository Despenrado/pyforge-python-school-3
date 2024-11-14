from sqlalchemy import update
from sqlalchemy.orm import class_mapper

from app.models.molecules import Molecule
from app.storage.postgres.db_connect import DBSessionLocal, get_db_session_local


class MoleculeStorage:
    def __init__(self, session_factory=DBSessionLocal):
        if session_factory is None:
            session_factory = get_db_session_local()
        self.session_factory = session_factory

    def save(self, molecule: Molecule):
        with self.session_factory() as session:
            print(molecule)
            session.add(molecule)
            session.commit()
            session.refresh(molecule)
        return molecule


    def get_all(self, limit: int = 1000, offset: int = 0):
        with self.session_factory() as session:
            query = session.query(Molecule).offset(offset).limit(limit)
            results = query.all()

        return results

    def get_by_id(self, molecule_id):
        with self.session_factory() as session:
            return session.query(Molecule).filter_by(id=molecule_id).first()


    def get_by_name(self, name):
        with self.session_factory() as session:
            return session.query(Molecule).filter_by(name=name).first()


    def get_by_smiles(self, smiles):
        with self.session_factory() as session:
            return session.query(Molecule).filter_by(smiles=smiles).first()


    def partial_update(self, molecule_id, molecule: Molecule):
        with self.session_factory() as session:
            columns = [column.key for column in class_mapper(Molecule).columns]
            update_data = {column: getattr(molecule, column) for column in columns}

            session.query(Molecule).filter(Molecule.id == molecule_id).update(
                update_data,
                synchronize_session="fetch"
            )
            session.commit()
            updated_molecule = session.query(Molecule).filter_by(id=molecule_id).first()
            return updated_molecule


    def delete(self, molecule_id):
        with self.session_factory() as session:
            molecule = session.query(Molecule).filter_by(id=molecule_id).first()
            if molecule:
                session.delete(molecule)
                session.commit()
            else:
                raise ValueError(f"Molecule with ID {molecule_id} not found.")
