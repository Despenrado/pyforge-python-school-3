from contextlib import asynccontextmanager

from fastapi import FastAPI

from app.routers import molecules
from app.storage.postgres.db_connect import db_init_engine, migrate
from app.storage.redis.redis_connect import connect_redis, disconnect_redis


@asynccontextmanager
async def lifespan(app: FastAPI):
    db_init_engine()
    migrate()
    connect_redis()
    yield

    disconnect_redis()

app = FastAPI(lifespan=lifespan)

app.include_router(molecules.mol_router, prefix="/api")