import os
from pathlib import Path
from matplotlib import font_manager
def configure_font(folder=None):
    folder=folder or os.environ.get('REPRO_FONT_DIR')
    if folder:
        folder=Path(folder)
        for name in ['arial.ttf','arialbd.ttf']:
            font_manager.fontManager.addfont(str(folder/name))
        return 'Arial'
    return 'DejaVu Sans'
