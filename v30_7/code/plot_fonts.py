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
    families={f.name for f in font_manager.fontManager.ttflist}
    return 'Liberation Sans' if 'Liberation Sans' in families else 'DejaVu Sans'
