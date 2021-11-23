# Disfuse-FORS-Analasis-Functions
Python Module Created for My Masters in Physics Project at Exeter University in Diffuse Fibre Optic Reflectence Spectral Analysis

ONLY FOR GOOGLE COLAB

Dependencies

numpy

scipy

scipy.interpolate

scipy.misc

scipy.signal

matplotlib.pyplot

glob

os

shutil

math

google.colab.auth

gspread

oauth2client.client.GoogleCredentials

Google Colab Headers to Import Function and Load Changes to Modual
```
### indicates different cell
!pip install GitPython
###
from google.colab import drive
drive.mount('/content/drive')
###
from git import Repo
import shutil
import sys
import importlib

#shutil.rmtree("/testing/")#uncomment for second run
Repo.clone_from("https://github.com/Haggi1181/Disfuse-FORS-Analasis-Functions.git", "/testing/")

sys.path.insert(1, "/testing/")

import MastersModual as mm
importlib.reload(mm)
###
```