#-to wipe out the container: docker image rm hemonc-spatial:latest
docker compose build --no-cache
docker compose up -d
#-open broswer to: http://localhost:8787
docker compose down
