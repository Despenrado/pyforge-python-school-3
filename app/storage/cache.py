import json

from app.storage.redis import redis_connect as redis_cache


def get_cached_result(key: str):
    result = redis_cache.get_cached_result(key)
    if result:
        return json.loads(result)
    return None

def set_cache(key: str, value, expiration: int = 60):
    redis_cache.set_cache(key,  json.dumps(value), expiration)