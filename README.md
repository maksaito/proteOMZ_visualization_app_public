# Application for Visualization of Marine Metaproteomics Data
This respository provides scripts and example data for the visualization app described in "Harnessing the power of Scientific Python to Investigate Metaproteomes and Biogeochemistry of the Central Pacific Ocean"
(Proceedings of the 17th Python in Science Conference, 2018) by Noelle Held, Jaclyn Saunders, Joe Futrelle, and Mak Saito. 

## Installing and running the application
The application can be executed by running `bokeh serve --show visualization_application.py`from within the repo. The application will open in your browser. Please note that the application is optimized for large format screens.

Note that the example data is scrambled and thus should not be used to make scientific conclusions. The complete dataset will be provided in an upcoming publication.

This application makes use of Bokeh (https://github.com/bokeh/bokeh) and pandas (https://github.com/pandas-dev/pandas). It is important to use the exact version of Bokeh in particular; see enviornment.yml.

## More Information
Description of Saito laboratory: http://www.whoi.edu/page.do?pid=36296

Data repository for proteOMZ research expedition: https://www.bco-dmo.org/project/685696

Description of proteOMZ research expedition: https://schmidtocean.org/cruise/investigating-life-without-oxygen-in-the-tropical-pacific/#team
