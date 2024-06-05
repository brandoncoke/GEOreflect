#/bin/sh
sudo docker stop $(sudo docker ps -a -q)
sudo docker rm $(sudo docker ps -a -q)
sudo  docker build -o - . --tag georeflect_shiny
sudo docker create georeflect_shiny
sudo docker container ls -a
sudo docker export --output="georeflect_shiny.tar.gz" georeflect_shiny
sudo docker container run --name georeflect_shiny --interactive --tty 

#sudo docker kill 
