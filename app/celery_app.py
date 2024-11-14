from celery import Celery

from app.storage.postgres.db_connect import db_init_engine
from app.storage.redis.redis_connect import connect_redis

celery = Celery(
    'worker',
    broker='redis://redis:6379/0',
    backend='redis://redis:6379/0',
    include=['app.services.tasks_molecules']
)

db_init_engine()
connect_redis()