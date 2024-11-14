import pytest
from app.models.molecules import Molecule
import app.storage.postgres.db_connect as db_connect
from app.storage.postgres.molecules import MoleculeStorage


@pytest.fixture(scope='module', autouse=True)
def setup_db():
    db_connect.db_init_engine('sqlite:///:memory:')
    db_connect.get_db_session_local()
    db_connect.migrate()
    yield

    db_connect.Base.metadata.drop_all(bind=db_connect.db_engine)

@pytest.fixture
def db_test_session():
    session = db_connect.DBSessionLocal()
    yield session
    session.close()

@pytest.fixture
def molecule_storage(db_test_session):
    storage = MoleculeStorage(session_factory=lambda: db_test_session)
    return storage


def test_save_molecule(molecule_storage):
    molecule = Molecule(name="Water", smiles="O")
    saved_molecule = molecule_storage.save(molecule)
    assert saved_molecule.id is not None
    assert saved_molecule.name == "Water"
    assert saved_molecule.smiles == "O"


def test_get_all_molecules(molecule_storage):
    molecules = molecule_storage.get_all()
    assert len(molecules) >= 1


def test_get_molecule_by_id(molecule_storage):
    molecule = Molecule(name="Ethanol", smiles="CCO")
    saved_molecule = molecule_storage.save(molecule)
    retrieved_molecule = molecule_storage.get_by_id(saved_molecule.id)
    assert retrieved_molecule.id == saved_molecule.id
    assert retrieved_molecule.name == "Ethanol"
    assert retrieved_molecule.smiles == "CCO"


def test_get_molecule_by_name(molecule_storage):
    molecule = Molecule(name="Methane", smiles="C")
    molecule_storage.save(molecule)
    retrieved_molecule = molecule_storage.get_by_name("Methane")
    assert retrieved_molecule.name == "Methane"
    assert retrieved_molecule.smiles == "C"


def test_get_molecule_by_smiles(molecule_storage):
    molecule = Molecule(name="Benzene", smiles="c1ccccc1")
    molecule_storage.save(molecule)
    retrieved_molecule = molecule_storage.get_by_smiles("c1ccccc1")
    assert retrieved_molecule.name == "Benzene"
    assert retrieved_molecule.smiles == "c1ccccc1"


def test_partial_update_molecule(molecule_storage):
    molecule = Molecule(name="Acetone", smiles="CC(=O)C")
    saved_molecule = molecule_storage.save(molecule)
    saved_molecule.name = "Propanone"
    updated_molecule = molecule_storage.partial_update(saved_molecule.id, saved_molecule)
    retrieved_molecule = molecule_storage.get_by_id(saved_molecule.id)
    assert retrieved_molecule.name == "Propanone"
    assert retrieved_molecule.smiles == "CC(=O)C"


def test_delete_molecule(molecule_storage):
    molecule = Molecule(name="Propane", smiles="CCC")
    saved_molecule = molecule_storage.save(molecule)
    molecule_storage.delete(saved_molecule.id)
    retrieved_molecule = molecule_storage.get_by_id(saved_molecule.id)
    assert retrieved_molecule is None
