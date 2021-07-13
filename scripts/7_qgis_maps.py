#!/usr/bin/env python3
# This scirpts read the QGIS project and export prepeared print layouts
# to 'plots' directory
import os
from qgis.core import (QgsProject, QgsLayoutExporter, QgsApplication)

QgsApplication.setPrefixPath("/usr", True)

gui_flag = False
app = QgsApplication([], gui_flag)

app.initQgis()

project_path = os.path.join('..', 'qgis','ruspolia.qgz')

project_instance = QgsProject.instance()
project_instance.setFileName(project_path)
project_instance.read()

manager = QgsProject.instance().layoutManager()
for i in manager.layouts():
    exporter = QgsLayoutExporter(i)
    exporter.exportToImage(os.path.join('..', 'plots',i.name()+'.png'),QgsLayoutExporter.ImageExportSettings())

app.exitQgis()
