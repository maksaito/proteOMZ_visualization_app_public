# Application for Visualization of Marine Metaproteomics Data
This respository provides scripts and example data for the visualization app described in "Harnessing the power of Scientific Python to Investigate Metaproteomes and Biogeochemistry of the Central Pacific Ocean"
(Proceedings of the 17th Python in Science Conference, 2018) by Noelle Held, Jaclyn Saunders, Joe Futrelle, and Mak Saito. 

## Installing and running the application
The application can be executed by running `bokeh serve --show visualization_application.py`from within the repo. The application will open in your browser. Please note that the application is optimized for large format screens.

Note that the example data is scrambled and thus should not be used to make scientific conclusions. The complete dataset will be provided in an upcoming publication.

This application makes use of Bokeh (https://github.com/bokeh/bokeh) and pandas (https://github.com/pandas-dev/pandas). It is important to use the exact version of Bokeh in particular; see enviornment.yml.

To make full use of the map plots, you will need to input your own Google API key into the visualization_application.py file (location is noted). To acquire a gmaps API key, see https://cloud.google.com/maps-platform/.

## More Information
Description of Saito laboratory: http://www.whoi.edu/page.do?pid=36296

Data repository for proteOMZ research expedition: https://www.bco-dmo.org/project/685696

Description of proteOMZ research expedition: https://schmidtocean.org/cruise/investigating-life-without-oxygen-in-the-tropical-pacific/#team

## Companion apps
Corresponding metaproteomics data: https://github.com/maksaito/proteOMZ_visualization_app_public

An application of datashader to the dataset: https://github.com/naheld/15000lines_datashader
