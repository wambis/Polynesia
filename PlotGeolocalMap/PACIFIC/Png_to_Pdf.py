#!/usr/bin/python
# -*- coding: utf-8 -*-
import obspy
import numpy as np
import os, sys
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter

from PIL import Image
from reportlab.pdfgen import canvas

#compute the percentiles for 2.5% and 97.5%





def png_to_pdf(png_path, pdf_path):
    image = Image.open(png_path)
    width, height = image.size

    c = canvas.Canvas(pdf_path, pagesize=(width, height))
    c.drawImage(png_path, 0, 0, width, height)
    c.save()

# Usage example
#png_path = "FigEvents_in_Chunk.png"
#pdf_path = "FigEvents_in_Chunk.pdf" 
#########################################
png_path = "Fig_Rapypath_OBS_in_Chunk.png"
pdf_path = "Fig_Rapypath_OBS_in_Chunk.pdf" 

#png_path = "FigStations_in_Chunk_L.png"
#pdf_path = "FigStations_in_Chunk_L.pdf" 

##############################
png_to_pdf(png_path, pdf_path) 
