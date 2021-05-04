import os
import numpy as np
import wx
from metaunidec.gui_elements.list_ctrls import ListCtrlPanel
from metaunidec.gui_elements.ud_cont_meta import main_controls
from metaunidec.gui_elements.ud_menu_meta import meta_menu
from unidec_modules.gui_elements import peaklistsort, mainwindow_base
import wx.lib.scrolledpanel as scrolled
from pubsub import pub
from unidec_modules import plot1d, plot2d



class MegaWindow(mainwindow_base.MainwindowBase):

    def __init__(self, parent, title, config=None, iconfile="logo.ico", *args, **kwargs):
        mainwindow_base.MainwindowBase.__init__(self, parent, title, config, iconfile)
        # wx.Frame.__init__(self, None, title=title)  # ,size=(200,-1))
        self.pres = parent
        if config is None:
            self.config = self.pres.eng.config
        else:
            self.config = config
            
        self.icon_path = iconfile