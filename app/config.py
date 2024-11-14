from dotenv import load_dotenv
import os

load_dotenv()

class Settings:
    DB_USER: str = os.getenv("DB_USER")
    DB_PASS: str = os.getenv("DB_PASS")
    DB_HOST: str = os.getenv("DB_HOST")
    DB_PORT: str = os.getenv("DB_PORT")
    DB_NAME: str = os.getenv("DB_NAME")
    INSTANCE_NAME: str = os.getenv("INSTANCE_NAME")

    @property
    def SQLALCHEMY_DATABASE_URL(self):
        db_url = os.getenv('SQLALCHEMY_DATABASE_URL')
        if db_url:
            return db_url
        else:
            return (
                f"postgresql://{self.DB_USER}:{self.DB_PASS}@{self.DB_HOST}:{self.DB_PORT}/{self.DB_NAME}"
            )


settings = Settings()