import json
import pytest
from fastapi.testclient import TestClient
from app.main import app
from app.models.molecules import Molecule
from app.storage.postgres import db_connect

from app.storage.postgres.molecules import MoleculeStorage
# import os
# os.environ['SQLALCHEMY_DATABASE_URL'] = 'sqlite:///:memory:'
# @asynccontextmanager
# async def lifespan_test(test_app: FastAPI):
#     from app.models.molecules import Molecule
#     db_connect.db_init_engine('sqlite:///:memory:')
#     db_connect.migrate()
#     yield
#
#     Base.metadata.drop_all(bind=db_connect.db_engine)
#
# test_app = FastAPI(lifespan=lifespan_test)
# test_app.include_router(molecules.mol_router, prefix="/api")
#
# client = TestClient(test_app)

client = TestClient(app)

@pytest.fixture(scope='module', autouse=True)
def setup_db():
    db_connect.db_init_engine('sqlite:///:memory:')
    db_connect.get_db_session_local()
    db_connect.migrate()
    yield

    db_connect.Base.metadata.drop_all(bind=db_connect.db_engine)

@pytest.fixture
def db_test_session(setup_db):
    session = db_connect.DBSessionLocal()
    yield session
    session.close()

@pytest.fixture
def molecule_storage(db_test_session):
    storage = MoleculeStorage(session_factory=lambda: db_test_session)
    return storage

def test_get_test_message():
    response = client.get("/api/")
    assert response.status_code == 200
    assert "instance_name" in response.json()

def test_create_molecule():
    molecule_data = {
        "name": "Water",
        "smiles": "O"
    }
    response = client.post("/api/molecules", json=molecule_data)
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == "Water"
    assert data["smiles"] == "O"
    assert "id" in data


