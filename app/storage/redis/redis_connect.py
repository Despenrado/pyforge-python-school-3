import logging
from curses import wrapper

import redis as r

logger = logging.getLogger(__name__)

redis_client = None


def connect_redis():
    logger.info("Redis connection established")
    global redis_client
    redis_client = r.Redis(host='redis', port=6379)
    try:
        redis_client.ping()
        logger.info("Redis connection established")
    except Exception as e:
        logger.error(f"Redis connection error: {e}")

def disconnect_redis():
    logger.info("Redis disconnecting...")
    redis_client.close()
    logger.info("Redis disconnected")


def get_cached_result(key: str):
    result = redis_client.get(key)
    return result

def set_cache(key: str, value: str, expiration: int = 60):
    redis_client.setex(key, expiration, value)