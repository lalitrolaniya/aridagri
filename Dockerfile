FROM r-base:4.4.0
RUN R -e "install.packages('aridagri', repos='https://cran.r-project.org')"
WORKDIR /home/analysis
CMD ["R"]
