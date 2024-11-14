FROM python:3.12-slim

WORKDIR /app

RUN mkdir app

COPY ../requirements.txt .
RUN pip install --no-cache-dir -r ./requirements.txt

COPY ../app /app/app
COPY ../utils /app/utils

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080"]