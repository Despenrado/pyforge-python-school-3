from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

from app.config import settings

Base = declarative_base()

db_engine = None
DBSessionLocal = None

def get_db_engine(database_url=None):
    if database_url is None:
        database_url = settings.SQLALCHEMY_DATABASE_URL
    return create_engine(database_url)

def db_init_engine(database_url=None):
    global db_engine
    db_engine = get_db_engine(database_url)

def get_db_session_local():
    global DBSessionLocal
    print('+++++++++++++++++++++++++++++++++DBSessionLocal++++++++++++++++++++++++++++++++++++++++++')
    print(DBSessionLocal)
    if DBSessionLocal is None:
        if db_engine is None:
            db_init_engine()
        DBSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=db_engine)
    print(DBSessionLocal)
    return DBSessionLocal

count = 0
def migrate():
    global count
    print(f'++++++++++++++++++++++++++++++++++db_engine {count}+++++++++++++++++++++++++++++++++++++++++')
    print(db_engine)

    count += 1
    if db_engine is None:
        db_init_engine()
    Base.metadata.create_all(bind=db_engine)
