

class Transaction: # сенеджер контекста

    def __enter__(self):
        print("Begin Transaction")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        print("End Transaction")
        return False


with Transaction() as t:
    raise Exception("Error")