
# syntax=docker/dockerfile:1
#sudo docker build .
#set image and working directory
FROM rocker/r-devel
WORKDIR /home/
COPY ../app_and_gpl570.R /home/
COPY ../example_RNAseq_DEG_frame.csv /home/
COPY ../percentile_matrix.RDS /home/
COPY ../percentile_matrix_p_value_RNAseq.RDS /home/

RUN R -e 'install.packages("shiny")'
RUN R -e 'install.packages("openxlsx")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("DT")'
RUN R -e 'install.packages("plotly")'
COPY app.R /app.R
EXPOSE 3838
CMD R -e 'shiny::runApp("app_and_gpl570.R", port = 3838, host = "0.0.0.0")'