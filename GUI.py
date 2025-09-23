#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021  Thomas Ferrari


import wx, wx.adv, wx.grid
import wx.lib.agw.advancedsplash as SPLASH
import subprocess, platform
import sys  # stdout
import csv
import threading
from pathlib import Path
from glob import glob

import chem
import shell
import output
from GUI_panels import *

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))  # set this file directory as CWD


# IMAGES
# img_QSARpy_splash = wx.Image('./pics/QSARpy_splash_black.png')
img_DQSAR_splash = wx.Image('./pics/D-QSAR_splash.png')
img_SARpy_logo = wx.Image('./pics/SARpy_logo.png')
img_snakes = wx.Image('./pics/SARpy_snakes.png')

# PATH
model_path = "../MODELS/"
dataset_path = "../DATASETS/"
database_path = "../DATABASES/"
export_path = "../EXPORTS/"
manual_html_path = "QSARpy2_Manual.html"
manual_pdf_path = "../QSARpy2_Manual.pdf"
DQSAR_path = 'DeltaQSAR_query.exe'


# frame_size = (1050, 700)
frame_size = (1280, 720)
line = '\u2500'
delta = '\u0394'


class SubThread(threading.Thread):
    def run(self):
        try:
            if self._target is not None:
                self.result = self._target(*self._args, **self._kwargs)
        finally:
            # Avoid a refcycle if the thread is running a function with
            # an argument that has a member that points to the thread.
            del self._target, self._args, self._kwargs


class StartTab(wx.Panel):

    def __init__(self, parent):
        super().__init__(parent)
        self.model = None  # DEBUG: maybe not needed....
        self.logo = wx.StaticBitmap(self, -1, img_snakes.ConvertToBitmap())
        self.st_welcome = wx.StaticText(self, -1, f"Welcome to {delta}-QSAR, ")
        self.st_toolkit = wx.StaticText(self, -1, "a knowledge discovery suite for automatic SAR/QSAR rules induction.")
        self.panel_SAR = StartPanel(self, " SARpy (classification) ")
        self.panel_QSAR = StartPanel(self, " QSARpy (regression) ")
        # self.st_rightarrow = wx.StaticText(self, -1, "\u2192")
        self.st_navigation = wx.StaticText(self, wx.ID_ANY, "Press [Tab] key to highlight navigation buttons")

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        font = wx.Font(14, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        self.st_welcome.SetFont(fonts.title)
        self.st_toolkit.SetFont(fonts.bigger)
        # self.st_rightarrow.SetFont(font)
        self.st_navigation.SetFont(fonts.italic)
        self.st_navigation.Hide()
        self.panel_SAR.button_load.SetFocus()

    def __do_layout(self):
        hbox_welcome = wx.BoxSizer(wx.HORIZONTAL)
        hbox_welcome.Add((10, -1))
        hbox_welcome.Add(self.st_welcome, flag=wx.ALIGN_BOTTOM | wx.BOTTOM, border=15)
        hbox_welcome.Add(self.st_toolkit, flag=wx.ALIGN_BOTTOM | wx.BOTTOM, border=15)
        vbox_logo = wx.BoxSizer(wx.VERTICAL)
        vbox_logo.Add(self.logo, flag=wx.ALIGN_RIGHT | wx.TOP | wx.RIGHT, border=10)
        hbox_welcome.Add(vbox_logo, 1)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.panel_SAR, 1, wx.RIGHT | wx.EXPAND, 20)
        # hbox.Add(self.st_rightarrow, 0, wx.ALIGN_CENTER_VERTICAL)
        hbox.Add(self.panel_QSAR, 1, wx.LEFT | wx.EXPAND, 20)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(hbox_welcome, 0, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=25)
        vbox.Add((-1, 10), 0)
        vbox.Add(hbox, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 40)
        vbox.Add((-1, 20))
        vbox.Add(self.st_navigation, flag=wx.LEFT, border=30)
        vbox.Add((-1, 20))
        vbox.Layout()
        self.SetSizer(vbox)

    def enableSAR(self):
        self.resetTab()
        self.panel_SAR.enable()
        self.panel_SAR.button.SetFocus()

    def enableQSAR(self):
        self.resetTab()
        self.panel_QSAR.enable()
        self.panel_QSAR.button.SetFocus()

    def resetTab(self):
        self.panel_SAR.enable(False)
        self.panel_QSAR.enable(False)
        self.panel_SAR.button_load.SetFocus()


class QSARpyTab(wx.Panel):
    console_columns = 90
    msgError = "Unexpected error!"
    base_labels = {'title': f"{delta}-QSAR",
                   'progress': "",
                   'left_button': "<<   Back",
                   'right_button': "Next   >>",
                   'path_label': "DATASET :",
                   'console_text': "",
                   'widget_title': "",
                   'link': "More...",
                   'button': "APPLY",
                   'checkbox': ""}

    def __init__(self, parent, custom_labels={}, title=None, text=None, index=None):
        super().__init__(parent)
        self.notebook = parent
        self.text = text
        self.index = index
        self.prev_page = None
        self.next_page = None

        self.logo = wx.StaticBitmap(self, -1, img_SARpy_logo.ConvertToBitmap())

        labels = QSARpyTab.base_labels.copy()
        labels.update(custom_labels)
        self.labels = labels
        self.title = wx.StaticText(self, label=title or labels['title'])
        self.progress = wx.StaticText(self, label=labels['progress'])

        cp = self.controlPanel = wx.Panel(self, style=wx.SUNKEN_BORDER)
        # DUMMY!
        # 1) create widget in the subclass (parent = self.controlPanel)
        # 2) use _setWidget in the subclass before _do_layout
        self.widget = wx.Panel(cp)
        self.link = wx.adv.HyperlinkCtrl(cp, label=labels['link'])
        self.button = wx.Button(cp, wx.ID_REFRESH, label=labels['button'])
        self.checkbox = wx.CheckBox(cp, label=labels['checkbox'])
        self.button.SetMinSize((200, 50))

        self.button_left = wx.Button(self, label=labels['left_button'])
        self.button_right = wx.Button(self, label=labels['right_button'])
        self.console = Console(self, labels)

        self.link.Bind(wx.adv.EVT_HYPERLINK, self.OnLink)
        self.button_right.Bind(wx.EVT_BUTTON, self.OnNext)
        self.button_left.Bind(wx.EVT_BUTTON, self.OnBack)

        self.__set_properties()

    def __set_properties(self):
        self.checkbox.Hide()
        self.checkbox.SetFont(fonts.smaller)
        self.title.SetFont(fonts.title)
        self.button_left.SetMinSize(button_minsize)
        self.button_right.SetMinSize(button_minsize)
        self.button.Disable()
        self.button_left.Disable()
        self.button_right.Disable()

    def _setWidget(self, panel):
        self.widget = panel

    def _do_layout(self):
        hbox_logo = wx.BoxSizer(wx.HORIZONTAL)
        hbox_logo.Add(self.logo)
        hbox_logo.Add(self.title, flag=wx.LEFT | wx.ALIGN_BOTTOM, border=10)
        hbox_logo.Add((1, -1), 1)
        hbox_logo.Add(self.progress, flag=wx.RIGHT | wx.ALIGN_BOTTOM, border=10)

        # CONTROL PANEL LAYOUT
        hbox_button = wx.BoxSizer(wx.HORIZONTAL)
        hbox_button.Add((1, -1), 1)
        hbox_button.Add(self.button, 2, flag=wx.EXPAND)
        hbox_button.Add((1, -1), 1)

        vbox_cp = wx.BoxSizer(wx.VERTICAL)
        vbox_cp.Add(self.widget, 2, flag=wx.EXPAND)
        # vbox_cp.Add((-1, 20))
        vbox_cp.Add((-1, 1), 1)  # expandable spacer
        vbox_cp.Add(self.link)
        vbox_cp.Add((-1, 1), 2)  # expandable spacer
        vbox_cp.Add(hbox_button, 1, flag=wx.BOTTOM | wx.EXPAND)
        vbox_cp.Add((-1, 1), 2)  # expandable spacer
        vbox_cp.Add(self.checkbox)
        vbox_cp.Add((-1, 5))
        hbox_cp = wx.BoxSizer(wx.HORIZONTAL)
        hbox_cp.Add((15, -1))  # padding left/right
        hbox_cp.Add(vbox_cp, 1, wx.EXPAND)
        hbox_cp.Add((15, -1))  # padding left/right
        hbox_cp.Layout()
        self.controlPanel.SetSizer(hbox_cp)

        # FOOTER BUTTONS
        hbox_footer = wx.BoxSizer(wx.HORIZONTAL)
        hbox_footer.Add(self.button_left, 1, wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        hbox_footer.Add(self.button_right, 1, wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

        vbox_left = wx.BoxSizer(wx.VERTICAL)
        vbox_left.Add(hbox_logo, flag=wx.ALL | wx.EXPAND, border=5)
        vbox_left.Add((-1, 10))
        vbox_left.Add(self.controlPanel, 1, flag=wx.EXPAND)
        vbox_left.Add(hbox_footer, flag=wx.TOP | wx.EXPAND, border=10)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(vbox_left, 2, flag=wx.ALL | wx.EXPAND, border=10)
        hbox.Add(self.console, 3, flag=wx.ALL | wx.EXPAND, border=10)
        hbox.Layout()
        self.SetSizer(hbox)

    def setTitle(self, text):
        self.title.SetLabel(text)

    def enableButton(self, enable=True):
        if enable:
            self.button.Enable()
            self.button.SetFocus()
        else:
            self.button.Disable()

    def resetTab(self, clear=True):
        self.title.SetLabel(self.labels['title'])
        self.enableButton(False)
        self.button_left.Disable()
        self.button_right.Disable()
        self.prev_page = self.next_page = None
        if clear:
            self.console.OnClear(prompt=self.labels['console_text'])

    def OnBack(self, event=None):
        self.notebook.ChangeSelection(self.notebook.GetSelection() - 1)
        # self.notebook.AdvanceSelection(forward=False)  # generates Page Changing event (it is vetoed)

    def OnNext(self, event=None):
        # self.notebook.AdvanceSelection()  # generates Page Changing event (it is vetoed)
        self.notebook.ChangeSelection(self.notebook.GetSelection() + 1)
        self.notebook.GetCurrentPage().button.SetFocus()

    def setOutput(self):  # INVECE DI METTERLO NEI BOTTONI, METTILO IN UN EVENTO CHANGE PAGE (c'era bug? anche se busy potevi cambiar pagina?)
        sys.stdout = self.console

    def writeDataset(self, filename):
        header = " DATASET: %s " % filename
        space = self.console_columns - len(header) - 1
        self.console.write("\n\n" + line + header + line * space + '\n', color=wx.BLUE)

    def writeRuleset(self, param_dict=None):
        self.console.write('\n\n' + line + " RULESET: ", color=wx.RED)
        self.console.write(str(self.model) + '\n')
        if param_dict is not None:
            self.console.write(f'  {param_dict}\n')

    def write(self, text):
        self.console.write('\n\n\n' + line + f' {text}\n', bold=True)

    # ONLY for Fragment and Predict Tabs
    def loadPath(self):
        with LoadDialog(self, 'dataset') as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                try:
                    with open(path) as f:
                        if path.suffix == '.csv':
                            headers = next(csv.reader(f))
                        elif path.suffix == '.sdf':
                            headers = []
                            line = f.readline()
                            while line and line[:4] != '$$$$':
                                if line.startswith('>'):
                                    start = line.find('<') + 1
                                    stop = line.find('>', start)
                                    headers.append(line[start:stop])
                                line = f.readline()
                        else:
                            raise OSError("Format '{}' is not supported (CSV and SDF only)".format(path.suffix))
                except Exception as e:
                    showMessage(dlg.msgErrorImport, exception=e, style=wx.ICON_ERROR)
                    raise
                else:
                    self.resetTab()
                    self.path = path
                    self.setupHeaders(headers)
                    self.console.setStatus(path.name)
                    self.writeDataset(path.name)
                    return True
            else:
                return False


class TabDialog(wx.Dialog):

    def __init__(self, parent, label, title=f'{delta}-QSAR settings'):
        super().__init__(parent, title=title)
        self.parent = parent
        self.title = wx.StaticText(self, label=label)
        self.title.SetFont(fonts.bold)
        self.widgets = {}
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def OnClose(self, event):
        self.Destroy()
        return wx.ID_CANCEL

    def _set_widgets(self, widget_dict):
        self.widgets = widget_dict

    def _do_layout(self, ratio=(1, 1)):  # set widgets before calling _do_layout
        button_sizer = self.CreateSeparatedButtonSizer(wx.OK)
        # wx.FindWindowById(wx.ID_OK, self).SetLabel('Apply')  # default apply/close doesn't work...
        # wx.FindWindowById(wx.ID_CANCEL, self).SetLabel('Close')
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.title, flag=wx.LEFT | wx.TOP, border=10)

        for key, widget in self.widgets.items():
            if key is None or widget is None:  # if no key or no widget, single line column
                vbox.Add(widget or wx.StaticText(self, label=key), 1, flag=wx.EXPAND | wx.ALL, border=20)
                continue
            label = wx.StaticText(self, label=key)
            hbox_widget = wx.BoxSizer(wx.HORIZONTAL)
            hbox_widget.Add(label, ratio[0], flag=wx.ALIGN_CENTER_VERTICAL)
            hbox_widget.Add(widget, ratio[1], flag=wx.ALIGN_CENTER_VERTICAL)
            vbox.Add(hbox_widget, 1, flag=wx.EXPAND | wx.ALL, border=20)
        vbox.Add(button_sizer, flag=wx.EXPAND)
        vbox.Add((-1, 10))

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(vbox, 1, flag=wx.EXPAND | wx.ALL, border=5)
        # sizer.Add((-1, 5))
        sizer.Layout()
        self.SetMinSize((400, 265))
        self.SetSizer(sizer)
        self.Fit()


class FragmentOptions(TabDialog):

    def __init__(self, parent):
        super().__init__(parent, "Fragmentation parameters:")
        self.spin = wx.SpinCtrl(self, value=str(parent.params['depth']), min=1, max=9)
        self.spin.SetMinSize(spinctrl_minsize)
        widgets = {"Fragmentation depth :": self.spin}
        if parent.type == 'QSAR':
            box = wx.StaticBox(self, label="QSAR parameters")
            sz = wx.StaticBoxSizer(box, wx.VERTICAL)
            self.QSAR = wx.CheckBox(box, label="Exclude outliers")
            self.QSAR.SetMinSize((-1, 50))
            sz.Add(self.QSAR, flag=wx.ALL, border=10)  # wx.SizerFlags().Border(wx.ALL, 30))
            self.QSAR.SetValue(parent.params['pop_outliers'])
            widgets[None] = sz
        self._set_widgets(widgets)
        self._do_layout()
        # self.SetSize((500, 500))  # DEBUG


# class BinarizeDialog(wx.Dialog):
#
#     checkbox_msg = "Detected {} classes!\n\n" \
#                    "Binarize the classes for optimal results: split the classes in two sets (checked/unchecked)" \
#                     " by ticking the relative checkbox."  #\nThen customize the labels."
#
#     textbox_msg = "Customize the binary labels:"
#
#     def __init__(self, parent, labels):
#         super().__init__(parent, title='Preprocessing data')
#         self.parent = parent
#         self.SetLayoutAdaptationMode(wx.DIALOG_ADAPTATION_MODE_ENABLED)  # AUTOMATIC SCROLL
#         checkbox_msg = wx.StaticText(self, label=self.checkbox_msg.format(len(labels)))
#         textbox_msg = wx.StaticText(self, label=self.textbox_msg)
#         button_sizer = self.CreateSeparatedButtonSizer(wx.OK)
#         self.checkbox_dict = {label: wx.CheckBox(self, label=label) for label in labels}
#         checked_label = wx.StaticText(self, label="checked  ")
#         unchecked_label = wx.StaticText(self, label="unchecked")
#         self.checked_textctrl = wx.TextCtrl(self, value="Positives")
#         self.unchecked_textctrl = wx.TextCtrl(self, value="Negatives")
#         self.ok = wx.FindWindowById(wx.ID_OK, self)
#         self.ok.Disable()
#
#         for cb in self.checkbox_dict.values():
#             cb.Bind(wx.EVT_CHECKBOX, self.OnChecking)
#
#         checkbox_msg.Wrap(500)
#
#         vbox = wx.BoxSizer(wx.VERTICAL)
#         vbox.Add(checkbox_msg)
#         vbox.Add((-1, 10))
#         for widget in self.checkbox_dict.values():
#             vbox.Add(widget, flag=wx.LEFT, border=5)
#
#         hbox_labels = wx.BoxSizer(wx.HORIZONTAL)
#         hbox_labels.Add(checked_label, 1)
#         hbox_labels.Add(self.checked_textctrl, 1)
#         hbox_textbox = wx.BoxSizer(wx.HORIZONTAL)
#         hbox_textbox.Add(unchecked_label, 1)
#         hbox_textbox.Add(self.unchecked_textctrl, 1)
#
#         vbox.Add((-1, 20))
#         vbox.Add(textbox_msg)
#         vbox.Add((-1, 10))
#         vbox.Add(hbox_labels, 1, wx.EXPAND | wx.LEFT, border=5)
#         vbox.Add(hbox_textbox, 1, wx.EXPAND | wx.LEFT, border=5)
#         vbox.Add((-1, 20))
#         vbox.Add(button_sizer, flag=wx.EXPAND)
#         sizer = wx.BoxSizer(wx.VERTICAL)
#         sizer.Add(vbox, 1, flag=wx.EXPAND | wx.ALL, border=10)
#
#         self.SetSizer(sizer)
#         self.Fit()
#         # self.SetSize((250, 500))
#
#     def OnChecking(self, event):
#         for cb in self.checkbox_dict.values():
#             if cb.GetValue():
#                 self.ok.Enable()
#                 return
#         self.ok.Disable()
#
#     def getDict(self):
#         checked_label, unchecked_label = self.checked_textctrl.GetValue(), self.unchecked_textctrl.GetValue()
#         result = {checked_label: set(), unchecked_label: set()}
#         for label, widget in self.checkbox_dict.items():
#             if widget.GetValue():
#                 result[checked_label].add(label)
#             else:
#                 result[unchecked_label].add(label)
#         return result


class FragmentTab(QSARpyTab):

    labels = {'title': f"{delta}-QSAR Extractor",
              'SAR_title': "SARpy Extractor",
              'QSAR_title': "QSARpy Extractor",
              'progress': "1/2",
              'path_label': "TRAINING SET :",
              'console_text': "\n Load a Training set from the File menu to start.",
              'link': "Advanced parameters...",
              'button': "FRAGMENT",
              'left_button': "Quick Export DB..."}

    SAR_PARAMETERS = {'depth': 4}
    QSAR_PARAMETERS = {'depth': 3, 'pop_outliers': False}

    def __init__(self, parent, text='   TRAINING   ', index=None):
        super().__init__(parent, custom_labels=FragmentTab.labels, text=text, index=index)
        self.trainingset = None
        self.model = None
        self.type = None
        self.params = None
        # CUSTOM WIDGET
        self.headerPanel = HeaderPanel(self.controlPanel, "SET TRAINING SET HEADERS:", outlier_checkbox=True)
        self._setWidget(self.headerPanel)
        # LAYOUT & PROPERTIES
        # self.button_right.Enable()  # DEBUG per vedere extract tab senza caricare trainingset
        # self.button_left.Hide()
        self.link.Enable(False)
        self._do_layout()
        # EVENTS
        self.headerPanel.button_load.Bind(wx.EVT_BUTTON, self.OnLoad)
        self.button.Bind(wx.EVT_BUTTON, self.OnFragment)

    def setType(self, type):
        self.type = type
        if type == 'SAR':
            self.title.SetLabel(self.labels['SAR_title'])
            self.params = self.SAR_PARAMETERS.copy()
            # self.button_left.SetLabel(self.labels['left_button'])
        elif type == 'QSAR':
            self.title.SetLabel(self.labels['QSAR_title'])
            self.params = self.QSAR_PARAMETERS.copy()
            # self.button_left.SetLabel('Quick Export DB...')
        else:
            raise

    def type(self):  # serve?
        if self.trainingset:
            return self.trainingset.type

    def resetTab(self, clear=True):
        # self.console.OnClear(prompt=FragmentTab.labels['console_text'])
        self.trainingset = None  # dovrebbe essere piu` efficiente e pulire anche self.model, self.type ecc.
        self.headerPanel.reset()
        self.link.Enable(False)
        super().resetTab(clear)

    def setupHeaders(self, headers):
        self.headerPanel.setHeaders(headers)
        self.headerPanel.setRequirements(ID=False, smiles=(self.path.suffix == '.csv'), target=True)

    # def _binarizeDialog(self, dataset):
    #     with BinarizeDialog(self, dataset.labels) as dlg:
    #         if dlg.ShowModal() == wx.ID_OK:
    #             dataset.binarize(dlg.getDict())

    def refresh(self):
        if self.trainingset:
            self.link.Enable()
            self.enableButton()
            self.button_right.Disable()
            self.button_left.Enable()  # "Quick export DB"
        else:
            self.link.Disable()
            self.enableButton(False)
            self.button_right.Disable()
            if self.type == 'QSAR':
                self.button_left.Disable()
        # if self.type == 'SAR':
        #     self.button_left.Enable()

    def OnLoad(self, event):
        wait = wx.BusyCursor()
        self.setOutput()
        idKey, smilesKey, targetKey = self.headerPanel.getHeaders()
        self.trainingset = None
        self.refresh()
        try:
            dataset = shell.loadDataset(self.path, targetKey, idKey, smilesKey)
            if self.type == 'QSAR':
                if not (trainingset := shell.makeQSARTrainingset(dataset)):
                    raise RuntimeWarning
            elif self.type == 'SAR':
                if not (trainingset := shell.makeSARTrainingset(dataset)):
                    raise RuntimeWarning
                labels = len(trainingset.labels)
                if labels > 2:
                    showMessage(f"Detected {labels} classes.\n"
                                f"\nEven if SARpy is multiclass, binary classification gives better results.",\
                                style=wx.ICON_INFORMATION)
                #     self._binarizeDialog(trainingset)
                #     print("\n   > CLASSES: ", trainingset)
            else:
                raise
        except ValueError:
            showMessage(f"'{targetKey}' NOT numeric", caption=f"{delta}-QSAR: Invalid operation", style=wx.ICON_EXCLAMATION)
        except RuntimeWarning:
            showMessage(f"'{smilesKey}' or '{targetKey}' values not valid",
                        caption=f"{delta}-QSAR: Invalid operation", style=wx.ICON_EXCLAMATION)
        except Exception as e:
            showMessage(self.msgError, exception=e, style=wx.ICON_ERROR)
            raise
        else:
            self.trainingset = trainingset
            self.refresh()
        finally:
            self.headerPanel.button_load.Disable()
            del wait

    def OnLink(self, event):
        with FragmentOptions(self) as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                # self.params['speed_up'] = dlg.checkbox.IsChecked()
                self.params['depth'] = dlg.spin.GetValue()
                if self.type == 'QSAR':  # QSAR
                    self.params['pop_outliers'] = dlg.QSAR.GetValue()
                # self.writeParameters(parameters)
                # if not self.trainingset or not self.trainingset.fragments or\
                #             self.depth != self.trainingset.fragmentation_level or \
                #             self.speed_up != (self.trainingset.max_atoms_fragmented is not None) or \
                #             (self.type() == 'SAR' and self.target != self.trainingset.target) or \
                #             (self.type() == 'QSAR' and self.pop_outliers == self.trainingset.excludeOutliers):
                #         self.writeParameters(parameters)
                #         # print("refactor needed")

    def OnFragment(self, event):
        self.setOutput()
        self.write(f"FRAGMENTATION:  Advanced parameters = {list(self.params.values())}")
        self.console.tc_logArea.Disable()
        wait = wx.BusyCursor()
        # if self.params.get('pop_outliers'):
        #     self.trainingset.popOutliers()
        speed_up = mainFrame.isSpeedUp()
        try:
            # kwargs = {'speed_up': speed_up, **self.params}
            # t = threading.Thread(target=shell.fragmentize, args=(self.trainingset,), kwargs=kwargs)
            # t.start()
            # t.join()
            shell.fragmentize(self.trainingset, speed_up=speed_up, **self.params)
        except Exception as e:
            self.trainingset.clear_fragmentation()
            showMessage(self.msgError, exception=e, style=wx.ICON_ERROR)
            raise
        else:
            if self.trainingset.fragments:
                self.button_right.Enable()
                self.button_right.SetFocus()
            else:
                self.button_right.Disable()
            self.next_page = None
        finally:
            del wait
            self.console.tc_logArea.Enable()

    def OnNext(self, event):
        index = self.notebook.GetSelection()
        if self.next_page is None:
            extractTab = ExtractTab(self.notebook, prev_page=self, text=self.text)
            self.next_page = extractTab
            self.notebook.InsertPage(index + 1, extractTab, extractTab.text)
            super().OnNext(event)
            self.notebook.RemovePage(index)
        else:
            self.notebook.InsertPage(index + 1, self.next_page, self.next_page.text)
            super().OnNext(event)
            self.notebook.RemovePage(index)

    def OnBack(self, event):  # If QSAR > EXPORT DATABASES
        if self.type == 'SAR':
            showMessage("This feature is not implemented yet in SARpy\n(QSARpy only!)", style=wx.ICON_INFORMATION)
            # super().OnBack(event)
        else:
            self.setOutput()
            self.console.tc_logArea.Disable()
            wait = wx.BusyCursor()
            try:
                shell.fragmentize(self.trainingset, depth=1)
                model = shell.makeQSARModel(self.trainingset, 0, 0, 0, verbose=False)
            except Exception as e:
                showMessage(self.msgError, exception=e, style=wx.ICON_ERROR)
            else:
                mainFrame.OnExportDB(model=model)
            finally:
                del wait
                self.console.tc_logArea.Enable()


class ExtractOptions(TabDialog):

    def __init__(self, parent):
        super().__init__(parent, "Extraction parameters:")
        self.spin_samples = wx.SpinCtrl(self, value=str(parent.params['min_samples']), min=2, max=99)
        self.spin_atoms = wx.SpinCtrl(self, value=str(parent.params['max_atoms']), min=1, max=20)
        self.spin_samples.SetMinSize(spinctrl_minsize)
        self.spin_atoms.SetMinSize(spinctrl_minsize)
        widgets = {"MIN sample size :": self.spin_samples,
                   "MAX atoms :": self.spin_atoms}

        box = wx.StaticBox(self)
        sz = wx.StaticBoxSizer(box, wx.VERTICAL)
        if parent.type() == 'SAR':
            box.SetLabel('SAR parameters')
            self.checkbox_H = wx.CheckBox(self, label='Explicit hydrogens')
            self.checkbox_exclusions = wx.CheckBox(self, label='Exclusions')
            self.checkbox_H.SetValue(parent.params['explicit_Hs'])
            self.checkbox_exclusions.SetValue(parent.params['exclusions'])
            # widgets['Explicit hydrogens :'] = self.checkbox_H
            # widgets['Exclusions :'] = self.checkbox_exclusions
            sz.Add(self.checkbox_H, flag=wx.ALL, border=10)
            sz.Add(self.checkbox_exclusions, flag=wx.ALL, border=10)
            widgets[None] = sz
        elif parent.type() == 'QSAR':
            box.SetLabel('QSAR parameters')
            self.checkbox_deep = wx.CheckBox(self, label='Deep learning')
            self.checkbox_deep.SetValue(parent.params['deep_learning'])
            sz.Add(self.checkbox_deep, flag=wx.ALL, border=10)
            widgets[None] = sz
            # widgets["Deep learning"] = self.checkbox_deep
        self._set_widgets(widgets)
        self._do_layout()


class ExtractTab(QSARpyTab):
    labels = {'title': f"{delta}-QSAR Extractor",
              'progress': "2/2",
              'path_label': "MODEL :",
              'console_text': "\n Calibrate the Model by adjusting the sliders.\n",
              'link': "Advanced parameters...",
              'button': "EXTRACT",
              'checkbox': " Enable CLUSTERING  (enhance model readability)",
              'right_button': f">>   send to {delta}-QSAR Predictor",
              'SAR_right_button': f">>   send to SARpy Predictor",
              'QSAR_right_button': f">>   send to QSARpy Predictor",
              }

    QSAR_widget_labels = {'title': "SET MAX ERROR:",
                          'left_label': "MAX accuracy",
                          'right_label': "MAX predictions",
                          # 'value': "|Value| :",
                          }

    SAR_PARAMETERS = {'min_samples': 4,
                      'max_atoms': 18,
                      'explicit_Hs': True,
                      'exclusions': False}
    QSAR_PARAMETERS = {'min_samples': 4,
                      'max_atoms': 7,
                      'deep_learning': False}

    def __init__(self, parent, prev_page, text=None):
        super().__init__(parent, custom_labels=ExtractTab.labels, title=prev_page.title.GetLabel(), text=text)
        self.path = None
        self.trainingset = prev_page.trainingset
        self.model = None
        if self.trainingset.type == 'SAR':
            self.button_right.SetLabel(self.labels['SAR_right_button'])
            self.widget = SAR_ExtractPanel(self.controlPanel, self.trainingset.labels)
            show_cluster_checkbox = False
        else:  # QSAR
            self.button_right.SetLabel(self.labels['QSAR_right_button'])
            self.widget = SliderPanel(self.controlPanel, ExtractTab.QSAR_widget_labels)
            self.widget.init(0, self.trainingset.abs_dev, self.trainingset.abs_dev * 2)
            show_cluster_checkbox = True
        # PARAMETERS
        if self.type() == 'SAR':
            self.params = self.SAR_PARAMETERS.copy()
        elif self.type() == 'QSAR':
            self.params = self.QSAR_PARAMETERS.copy()
        # HIERARCHY
        self.prev_page = prev_page
        # CUSTOM WIDGET
        self._setWidget(self.widget)
        # LAYOUT & PROPERTIES
        self.checkbox.Show(show=show_cluster_checkbox)
        self.checkbox.Disable()  # CHKBX
        self.button.Enable()
        self.button_left.Enable()
        self.button_right.SetId(wx.ID_FORWARD)
        self._do_layout()
        # EVENTS
        self.button.Bind(wx.EVT_BUTTON, self.OnExtract)
        self.checkbox.Bind(wx.EVT_CHECKBOX, self.OnFlowChanged)

    def type(self):
        return self.trainingset.type

    def OnLink(self, event):
        with ExtractOptions(self) as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                self.params['min_samples'] = dlg.spin_samples.GetValue()
                self.params['max_atoms'] = dlg.spin_atoms.GetValue()
                if self.type() == 'QSAR':
                    self.params['deep_learning'] = dlg.checkbox_deep.IsChecked()
                elif self.type() == 'SAR':
                    self.params['explicit_Hs'] = dlg.checkbox_H.IsChecked()
                    self.params['exclusions'] = dlg.checkbox_exclusions.IsChecked()

    def OnExtract(self, event):
        wait = wx.BusyCursor()
        self.setOutput()

        # # SAR: ANALYZE DECOMPOSITION (if 'EXTRACT" again without changing any param)
        # if self.type() == 'SAR' and self.model and self.model.fragments is None and\
        #         self.params.items() | self.widget.getDict().items() <= self.model.items():
        #     try:
        #         shell.analyzeDecomposition(self.model)
        #         self.writeRuleset()
        #         self.console.setStatus(str(self.model))
        #     except Exception as e:
        #         showMessage(self.msgError, exception=e, style=wx.ICON_ERROR)
        #         raise
        #     finally:
        #         del wait
        #     return

        self.write(f"EXTRACTION:  Advanced parameters = {list(self.params.values())}")
        try:
            if self.type() == 'SAR':
                widget_params = self.widget.getDict()
                model = shell.makeSARModel(self.trainingset, **self.params, **widget_params)
                shell.analyzeDecomposition(model)
                # model_constructor = shell.makeSARModel
            else:  # QSAR
                widget_params = {'max_error': self.widget.getValue()}
                model = shell.makeQSARModel(self.trainingset, **self.params, **widget_params)
                # model_constructor = shell.makeQSARModel
            # thread = SubThread(target=model_constructor, args=(self.trainingset,), kwargs=self.params)
            # thread.start()
            # thread.join()
            # model = thread.result
            # model = model_constructor(self.trainingset, **self.params)
            if model.ruleset:
                shell.validateOnTrain(model)
        except Exception as e:
            showMessage(self.msgError, exception=e, style=wx.ICON_ERROR)
            raise
        else:
            self.model = model
            self.console.setStatus(str(self.model))
            self.writeRuleset({k: widget_params[k] for k in ('max_error', 'beta') if k in widget_params})
            # if self.type() == 'SAR' and len(model.trainingset.labels) == 2:  # BINARY DATASET
            #     dialog = "BINARY CLASSIFICATION: hit [EXTRACT] again to analyze the structural decomposition"
            #     if self.model.isMonoTarget():
            #         dialog += "\n  (increase sensitivity)"
            #     self.write(dialog)
            if model.ruleset:
                self.checkbox.Enable()  # CHKBX
                self.button_right.Enable()
                self.button_right.SetFocus()
            else:
                self.checkbox.Disable()  # CHKBX
                self.button_right.Disable()
            self.next_page = None
            event.Skip()
        finally:
            del wait

    def OnBack(self, event):
        index = self.notebook.GetSelection()
        self.notebook.InsertPage(index, self.prev_page, self.prev_page.text)
        super().OnBack(event)
        self.notebook.RemovePage(index + 1)

    def OnNext(self, event):
        index = self.notebook.GetSelection()
        if self.checkbox.IsChecked():
            if self.next_page is None:
                generalizeTab = GeneralizeTab(self.notebook, prev_page=self, text=self.text)
                self.next_page = generalizeTab
                self.notebook.InsertPage(index + 1, generalizeTab, generalizeTab.text)
                super().OnNext(event)
                self.notebook.RemovePage(index)
            else:
                self.notebook.InsertPage(index + 1, self.next_page, self.next_page.text)
                super().OnNext(event)
                self.notebook.RemovePage(index)
        else:
            next_page = self.notebook.GetPage(index + 1)
            next_page.importModel(self.model)
            super().OnNext(event)
            next_page.showLoadDatasetDialog()

    def OnFlowChanged(self, event=None):  # QSARpy only
        if self.checkbox.IsChecked():
            self.progress.SetLabel('2/3')
            self.button_right.SetLabel("Clustering   >>")
        else:
            self.progress.SetLabel(ExtractTab.labels['progress'])
            self.button_right.SetLabel(ExtractTab.labels['QSAR_right_button'])
        self.button_right.SetFocus()


class GeneralizeTab(QSARpyTab):
    labels = {'title': "QSARpy Extractor",
              'progress': "3/3",
              'path_label': "MODEL :",
              'console_text': "\n Rationalize the Model by clustering the modulators. ",
              'link': "Automatic setting >",
              'button': "CLUSTER MODEL",
              'right_button': f">>   send to QSARpy Predictor"}

    widget_labels = {'title': "SET MODEL TOLERANCE:",
                     'left_label': "STRICT",
                     'left_label_2': "",
                     'right_label': "FEWER and",
                     'right_label_2': "BROADER rules",
                     'value': "Value : \u00B1",
                     'text': "(account minor endpoint variations as experimental error)"}

    def __init__(self, parent, prev_page, text=None):
        super().__init__(parent, custom_labels=GeneralizeTab.labels, text=text)
        self.original_model = prev_page.model  # unclustered model
        self.model = prev_page.model
        # HIERARCHY
        self.prev_page = prev_page
        # CUSTOM WIDGET
        self.widget = SliderPanel(self.controlPanel, GeneralizeTab.widget_labels)
        self.widget.init(0, self.model.max_error / 4, self.model.max_error / 2)
        self._setWidget(self.widget)
        # PROPERTIES & LAYOUT
        self.console.setStatus(str(self.model))
        self.link.Hide()
        self.button.Enable()
        self.button_left.Enable()
        self.button_right.SetId(wx.ID_FORWARD)
        self._do_layout()
        # EVENTS
        self.button.Bind(wx.EVT_BUTTON, self.OnGeneralize)

    def OnLink(self, event):
        pass

    def OnGeneralize(self, event):
        wait = wx.BusyCursor()
        self.setOutput()
        tolerance = self.widget.getValue()
        try:
            clusterModel = shell.getClustered(self.original_model, tolerance)
            if clusterModel.ruleset:
                shell.validateOnTrain(clusterModel)
        except Exception as e:
            showMessage(self.msgError, exception=e, style=wx.ICON_ERROR)
            raise
        else:
            self.model = clusterModel
            self.console.setStatus(str(self.model))
            self.writeRuleset({'tolerance': tolerance})
            if self.model.ruleset:
                self.button_right.Enable()
                self.button_right.SetFocus()
            else:
                self.button_right.Disable()
            event.Skip()
        finally:
            del wait

    def OnBack(self, event):
        index = self.notebook.GetSelection()
        self.notebook.InsertPage(index, self.prev_page, self.prev_page.text)
        super().OnBack(event)
        self.notebook.RemovePage(index + 1)

    def OnNext(self, event):
        next_page = self.notebook.GetPage(self.notebook.GetSelection() + 1)
        next_page.importModel(self.model)
        super().OnNext(event)
        next_page.showLoadDatasetDialog()

    def OnRefreshLabels(self):
        self.button_right.SetLabel(GeneralizeTab.labels['right_button'])


class PredictTab(QSARpyTab):
    labels = {'title': f"{delta}-QSAR Predictor",
              'SAR_title': "SARpy Predictor",
              'QSAR_title': "QSARpy Predictor",
              'right_button': "PREVIEW Predictions",
              'path_label': "TEST SET :",
              'console_text': "\n Load a Test set from the File menu to start.",
              'link': "Open model...",
              'button': "PREDICT"}

    def __init__(self, parent, text="   TESTING   ", index=None):
        super().__init__(parent, custom_labels=PredictTab.labels, text=text, index=index)
        self.path = None
        self.testset = None
        self.model = None
        # CUSTOM WIDGET
        self.headerPanel = HeaderPanel(self.controlPanel, "SET TEST SET HEADERS:")
        self._setWidget(self.headerPanel)
        # PROPERTIES & LAYOUT
        self.link.SetId(wx.ID_REFRESH)
        self._do_layout()
        # EVENTS
        self.headerPanel.button_load.Bind(wx.EVT_BUTTON, self.OnLoad)
        self.button.Bind(wx.EVT_BUTTON, self.OnPredict)
        self.button_right.Bind(wx.EVT_BUTTON, self.OnViewPredictions)

    def resetTab(self, clear=False):
        # self.console.OnClear(prompt=PredictTab.labels['console_text'])
        self.testset = None  # dovrebbe essere piu` efficiente e pulire anche self.model, self.path ecc.
        self.headerPanel.reset()
        super().resetTab(clear)

    def showLoadDatasetDialog(self):
        if self.path is None:
            if showMessage("Load a Test set?", style=wx.OK | wx.CANCEL) == wx.ID_OK:
                self.loadPath()

    def importModel(self, model, from_file=False):
        if model is not self.model:
            # if model['type'] != self.model['type']:
                
            if model['type'] == 'SAR':
                self.title.SetLabel(self.labels['SAR_title'])
            if model['type'] == 'QSAR':
                self.title.SetLabel(self.labels['QSAR_title'])
            self.model = model
            if not from_file:
                self.button_left.Enable()
                self.writeRuleset()
            super().writeRuleset()
            if self.testset:
                self.testset.reset()  # clear predictions
                self.enableButton()

    def writeRuleset(self, filename=None, parameters=None):
        model = parameters or self.model
        # model = self.model
        header = f" {model['type']} MODEL: "
        if filename is not None:
            header += filename + ' '
        space = self.console_columns - len(header) - 1
        self.console.write("\n\n" + line + header + line * space + '\n', color=wx.RED)
        self.console.write("   TRAINING SET:  < %s >\n" % model['train_name'])
        frag_parameters = [model['depth']]
        extr_parameters = [model['min_samples'], model['max_atoms']]
        if model['type'] == 'SAR':
            extr_parameters.append(model['explicit_Hs'])
            extr_parameters.append(model['exclusions'])
        else:
            frag_parameters.append(model['pop_outliers'])
            extr_parameters.append(model['deep_learning'])
        options = {'Chemical library': model['chem_library'], 'Speed-up': model['speed_up']}
        self.console.write(f"   Options = {options}\n")
        self.console.write(f"   Fragmentation parameters = {frag_parameters}\n")
        self.console.write(f"   Extraction parameters = {extr_parameters}\n")
        if model['type'] == 'SAR':
            sliders = [model['max_error'], model['beta']]
            self.console.write(f"   Target = {model['mono_target']}\n")
            self.console.write(f"   Sliders = {sliders}")
        else:
            self.console.write(f"   MAX error = {model['max_error']}\n")
            if model.get('tolerance') is not None:
                self.console.write(f"   Tolerance = {model['tolerance']}\n")
        # self.console.write("\n   RULESET: ", color=wx.RED)
        # self.console.write(str(self.model) + '\n')

    def setupHeaders(self, headers):
        self.headerPanel.setHeaders(headers)
        self.headerPanel.setRequirements(ID=False, smiles=(self.path.suffix == '.csv'), target=False)

    def OnLoad(self, event=None):
        wait = wx.BusyCursor()
        self.setOutput()
        idKey, smilesKey, targetKey = self.headerPanel.getHeaders()
        try:
            dataset = shell.loadDataset(self.path, targetKey, idKey, smilesKey)
            testset = shell.Testset(dataset)
        except Exception as e:
            showMessage(self.msgError, exception=e, style=wx.ICON_ERROR)
            raise
        else:
            self.testset = testset
        finally:
            self.headerPanel.button_load.Disable()
            if self.model:
                self.enableButton()
            else:
                self.link.SetFocus()
            del wait

    def OnLink(self, event):
        with LoadDialog(self, 'ruleset') as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                wx.BeginBusyCursor()
                self.setOutput()
                try:
                    model_stamp = shell.loadModel(path)
                    parameters = model_stamp['PARAMETERS']
                    mainFrame.updateOptions(parameters)
                    self.writeRuleset(path.name, parameters=parameters)
                    self.console.write("\n\n\n RE-EXTRACTING MODEL...", color=wx.RED)
                    model = shell.re_extractModel(model_stamp)
                except Exception as e:
                    showMessage(dlg.msgErrorImport, exception=e, style=wx.ICON_ERROR)
                    raise
                else:
                    if model is None:
                        self.console.write("\n ERROR!", color=wx.RED)
                        if self.model is not None:
                            self.console.write("  - Restoring previous model")
                            # self.writeImportedRuleset(current_model, "*RESTORED MODEL*")
                        else:
                            self.console.write("  - Model can't be loaded")
                    else:
                        super().resetTab(clear=False)  # disable buttons
                        self.importModel(model)
                        self.showLoadDatasetDialog()
                        event.Skip()
                finally:
                    wx.EndBusyCursor()

    def OnPredict(self, event):
        wait = wx.BusyCursor()
        self.setOutput()
        mode = self.model.trainingset.type
        # self.testset.reset()  # QSAR: do not reset; SAR: button disabled
        try:
            if self.testset.predictions:  # incremental prediction if already predicted
                if self.testset.unpredicted:
                    shell.predict_further(self.testset, self.model)
                    predicted_further = True
            else:
                predicted_further = False
                if mode == 'QSAR':
                    shell.fragmentize(self.testset, depth=self.model.depth)  # Here I fragment also potential identities...
                shell.predict(self.testset, self.model)
        except Exception as e:
            showMessage(self.msgError, exception=e, style=wx.ICON_ERROR)
            raise
        else:
            if self.testset.predictions:
                if self.headerPanel.getHeaders(2):  # IF NOT blind test set
                    try:
                        self.model.printStatistics(self.testset)
                    except ValueError:
                        showMessage("Prediction statistics are not shown because the selected endpoint is not numeric", style=wx.ICON_ERROR)
                self.button_left.Enable()
                self.button_right.Enable()
                self.button_right.SetFocus()
                if self.testset.unpredicted:  # predict further
                    if mode == 'SAR':
                        if not predicted_further and self.model.fragments:
                            dialog = "BINARY CLASSIFICATION: hit [PREDICT] again to predict further using QSARpy decomposition"
                            if self.model.isMonoTarget():
                                dialog += "\n  (increase sensitivity)"
                            self.write(dialog)
                        else:
                            self.enableButton(False)
                    elif mode == 'QSAR':
                        self.write(f"{len(self.testset.unpredicted)} unpredicted: hit [PREDICT] to predict further")
                else:
                    self.enableButton(False)
            else:
                self.button_right.Disable()
                self.enableButton(False)  # if unpredicted stop, don't predict further, messes uo the workflow...  I SHOULD FIX THIS!!!!!!!!!!!!!!!!!!!!!
            # if self.testset.unpredicted:  # predict further
            #     self.write("{} unpredicted: hit [PREDICT] to predict further".format(len(self.testset.unpredicted)))
            # else:
            #     self.enableButton(False)
            event.Skip()
        finally:
            del wait

    def OnNext(self, event):
        event.Skip()

    def OnViewPredictions(self, event):
        dataset = self.testset
        table = [['SMILES', 'prediction']]
        for smiles, pred_dict in dataset.predictions.items():
            table.append([smiles, pred_dict['prediction']])
        title = f"{delta}-QSAR predictions on dataset: {dataset.path}"
        with GridDialog(self, table, title) as dlg:
            dlg.ShowModal()


class GridDialog(wx.Dialog):

    def __init__(self, parent, table, title):
        super().__init__(parent, title=title, style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE)
        grid = wx.grid.Grid(self, 0)
        headers = table.pop(0)
        grid.CreateGrid(len(table), len(headers))
        for i, label in enumerate(headers):
            grid.SetColLabelValue(i, label)
        for row in range(len(table)):
            for col in range(len(headers)):
                grid.SetCellValue(row, col, str(table[row][col]))
        # grid.AutoSizeColumns(False)
        size = parent.GetSize()[1]
        grid.SetColSize(0, size - 200)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(grid, 1, wx.EXPAND)
        sizer.Layout()
        self.SetSizer(sizer)
        self.SetMinSize((size, size))
        #self.SetMaxSize(frame_size)
        #self.SetInitialSize()


class LoadDialog(wx.FileDialog):
    msgErrorImport = "\nAn error has been encountered while loading the selected file."
    msgOkImport = "\nFile '%s' has been loaded."
    keys = {'dataset': ("LOAD Dataset", dataset_path, "Comma Separated Values (*.csv)|*.csv|Structure Data File (*.sdf)|*.sdf"),
            'ruleset': ("LOAD Model", model_path, "JSON file (*.json)|*.json")}
    last_path_dict = {'dataset': None,
                      'ruleset': None
                      }

    def __init__(self, parent, key):
        try:
            defaultDir = Path(LoadDialog.last_path_dict[key])
        except:
            defaultDir = Path(self.keys[key][1])
            defaultDir.mkdir(exist_ok=True)
        super().__init__(parent, self.keys[key][0], str(defaultDir), wildcard=self.keys[key][2], style=wx.FD_OPEN)
        self.key = key

    def GetPath(self):
        path = Path(super().GetPath())
        LoadDialog.last_path_dict[self.key] = str(path.parent)  # save last directory
        return path


class SaveDialog(wx.FileDialog):
    msgErrorExport = "\nAn error has been encountered while saving the selected file."
    msgOkExport = "\nFile '%s' has been saved."
    keys = {'CSV': ("SAVE", export_path, "Comma Separated Values (*.csv)|*.csv", '.csv'),
            'ruleset': ("SAVE Model", model_path, "JSON file (*.json)|*.json", '.json'),
            'text': ("SAVE log", export_path, "Text file (*.txt)|*.txt", '.txt'),
            'database': ("EXPORT Database", database_path, "SQLite database (*.db)|*.db", '.db')}

    def __init__(self, parent, key):
        defaultDir = Path(self.keys[key][1])
        defaultDir.mkdir(exist_ok=True)
        super().__init__(parent, self.keys[key][0], str(defaultDir), wildcard=self.keys[key][2], style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        self.ext = self.keys[key][3]

    def GetPath(self):  # add suffix if not present
        path = super().GetPath()
        if path[-len(self.ext):] != self.ext:
            path += self.ext
        return Path(path)


class DBDialog(TabDialog):

    labels = (("Name :", 'Example: "Water solubility"'),
              ("Unit :", 'Example: "mg/l" or blank'),
              (None,  """Choose how to map the "endpoint shift"
                      \n Cyclic colorbar (DEFAULT): 
                      \n\t\u2212 [ BLU / RED ] +
                      \n Diverging colormap (e.g., TOXIC ENDPOINTS):
                      \n\t\u2212 [ RED / GREEN ] +\t\t( positive shift \u2192 "safe" ) 
                      \n\t\u2212 [ GREEN / RED ] +\t\t( negative shift \u2192 "safe" ) \n""".expandtabs(4)),
              )

    def __init__(self, parent):
        super().__init__(parent, "ENDPOINT information:", 'Exporting DATABASE...')
        self.name = wx.TextCtrl(self)
        self.unit = wx.TextCtrl(self)
        size = (250, -1)
        self.name.SetMinSize(size)
        self.unit.SetMinSize(size)

        cb_labels = ('twilight_shifted', 'RdYlGn', 'RdYlGn_r')
        cb_images = (wx.Image('./pics/colorbar_cyclic.png', wx.BITMAP_TYPE_ANY),
                     wx.Image('./pics/colorbar.png', wx.BITMAP_TYPE_ANY),
                     wx.Image('./pics/colorbar_reversed.png', wx.BITMAP_TYPE_ANY))

        self.rb_dict = {}
        self.cb_img_dict = {}
        value = True
        for label, img in zip(cb_labels, cb_images):
            self.rb_dict[label] = wx.RadioButton(self)
            self.rb_dict[label].SetValue(value)
            value = False
            scaled_img = img.Scale(size[0] + 20, 20)
            self.cb_img_dict[label] = wx.StaticBitmap(self, -1, wx.Bitmap(scaled_img))

        # COLORMAP
        vbox = wx.BoxSizer(wx.VERTICAL)
        for key in cb_labels:
            hbox = wx.BoxSizer(wx.HORIZONTAL)
            hbox.Add(self.rb_dict[key])
            hbox.Add(wx.StaticText(self, label="\u2212 "))
            hbox.Add(self.cb_img_dict[key])
            hbox.Add(wx.StaticText(self, label=" \u002B"))
            vbox.Add(hbox)
        sz = wx.StaticBoxSizer(wx.VERTICAL, self, "Colormap")
        sz.Add(vbox, flag=wx.ALL, border=10)

        widgets = {}
        for labels, widget in zip(self.labels, (self.name, self.unit, sz)):
            hbox = wx.BoxSizer(wx.HORIZONTAL)
            hbox.Add(widget)
            try:
                widget.SetHint(labels[1])
            except AttributeError:  # ADD TOOLTIP
                wx.ToolTip.SetAutoPop(20000)  # timer (ms)
                img_info = wx.StaticBitmap(self, bitmap=wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, size=(20, 20)))
                img_info.SetToolTip(labels[1])
                hbox.Add(img_info, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
            widgets[labels[0]] = hbox

        self._set_widgets(widgets)
        self._do_layout((1, 3))
        self.SetMinSize((-1, 400))

    def getdata(self):
        cb = None
        for key, rb in self.rb_dict.items():
            if rb.GetValue():
                cb = key
                break
        return {'endpoint': self.name.GetValue(),
                'unit': self.unit.GetValue(),
                'colorbar': cb}


class QSARpyMenuBar(wx.MenuBar):

    def __init__(self):
        super().__init__()
        # File Menu
        file_menu = wx.Menu()
        file_menu.Append(1, "Load SAR Training set...")
        file_menu.Append(0, "Load QSAR Training set...")
        file_menu.Append(2, "Load Test set...")
        file_menu.AppendSeparator()
        file_menu.Append(wx.ID_OPEN, "Open Model...")
        file_menu.AppendSeparator()
        file_menu.Append(wx.ID_EXIT, "Quit")
        self.Append(file_menu, '&File')

        # QSARpy Extractor Save Menu (train)
        save_menu = wx.Menu()
        save_menu.Append(wx.ID_CONTEXT_HELP, "Context: Extractor")
        save_menu.Enable(wx.ID_CONTEXT_HELP, False)
        save_menu.AppendSeparator()
        save_menu.AppendSeparator()
        save_menu.Append(wx.ID_SAVE, "Save Model...")
        save_menu.Append(5, "Export Predictions (Train)...")
        save_menu.AppendSeparator()
        save_menu.Append(7, "Export Model sheet...")
        save_menu.Append(9, "Export Clustered Model sheet...")
        save_menu.Append(8, "Export endpoint DATABASE...")
        self.TR_save_menu = save_menu
        self.Append(self.TR_save_menu, '&Save')

        # SARpy Extractor Save Menu (train)
        save_menu = wx.Menu()
        save_menu.Append(wx.ID_CONTEXT_HELP, "Context: Extractor")
        save_menu.Enable(wx.ID_CONTEXT_HELP, False)
        save_menu.AppendSeparator()
        save_menu.AppendSeparator()
        save_menu.Append(wx.ID_SAVE, "Save Model...")
        save_menu.Append(7, "Export Model sheet...")
        save_menu.Append(5, "Export Predictions (Train)...")
        self.SAR_TR_save_menu = save_menu
        # self.Append(self.SAR_TR_save_menu, '&Save')

        # QSARpy Predictor Save Menu (test)
        menu = wx.Menu()
        menu.Append(wx.ID_CONTEXT_HELP, "Context: Predictor")
        menu.Enable(wx.ID_CONTEXT_HELP, False)
        menu.AppendSeparator()
        menu.AppendSeparator()
        menu.Append(wx.ID_SAVE, "Save Model...")
        menu.Append(6, "Save Predictions (Test)...")
        menu.AppendSeparator()
        menu.Append(7, "Export Model sheet...")
        menu.Append(9, "Export Clustered Model sheet...")
        menu.Append(8, "Export endpoint DATABASE...")
        self.TE_save_menu = menu

        # SARpy Predictor Save Menu (test)
        menu = wx.Menu()
        menu.Append(wx.ID_CONTEXT_HELP, "Context: Predictor")
        menu.Enable(wx.ID_CONTEXT_HELP, False)
        menu.AppendSeparator()
        menu.AppendSeparator()
        menu.Append(wx.ID_SAVE, "Save Model...")
        menu.Append(7, "Export Model sheet...")
        menu.Append(6, "Save Predictions (Test)...")
        self.SAR_TE_save_menu = menu

        # Options Menu
        opt_menu = wx.Menu()
        ob = wx.MenuItem(opt_menu, id=10, text='Open Babel', kind=wx.ITEM_RADIO)
        rdkit = wx.MenuItem(opt_menu, id=11, text='RDKit', kind=wx.ITEM_RADIO)
        opt_menu.Append(ob)
        opt_menu.Append(rdkit)
        opt_menu.AppendSeparator()
        speed_up = wx.MenuItem(opt_menu, id=12, text='Speed-up', kind=wx.ITEM_CHECK)
        opt_menu.Append(speed_up)
        speed_up.Check()
        self.Append(opt_menu, '&Options')

        # Database Menu
        view_menu = wx.Menu()
        view_menu.Append(30, "Query database (QSAR only)...")
        self.delta_menu = view_menu
        self.Append(self.delta_menu, "Database")

        # Help Menu
        help_menu = wx.Menu()
        help_menu.Append(22, "SARpy Manual")
        help_menu.Append(20, "QSARpy Manual")
        help_menu.Append(21, f"About")
        self.Append(help_menu, "&Help")

        # QSARpy PREDICTOR Save Menu

    def display(self, context, mode):
        menu_idx = self.FindMenu('Save')
        if mode == 'SAR':
            if context == 'train' and self.GetMenu != self.SAR_TR_save_menu:
                self.Remove(menu_idx)
                self.Insert(menu_idx, self.SAR_TR_save_menu, '&Save')
            elif context == 'test' and self.GetMenu != self.SAR_TE_save_menu:
                self.Remove(menu_idx)
                self.Insert(menu_idx, self.SAR_TE_save_menu, '&Save')
        elif mode == 'QSAR':
            if context == 'train' and self.GetMenu != self.TR_save_menu:
                self.Remove(menu_idx)
                self.Insert(menu_idx, self.TR_save_menu, '&Save')
            elif context == 'test' and self.GetMenu != self.TE_save_menu:
                self.Remove(menu_idx)
                self.Insert(menu_idx, self.TE_save_menu, '&Save')

    def enableSave(self, enable=True):
        menu_idx = self.FindMenu('Save')
        self.EnableTop(menu_idx, enable)


class QSARpyFrame(wx.Frame):

    def __init__(self, parent, title):
        super().__init__(parent, title=title, size=frame_size)
        self.SetIcon(wx.Icon('./pics/QSARpy.ico', wx.BITMAP_TYPE_ICO))

        # MENU
        self.SetMenuBar(QSARpyMenuBar())
        self.MenuBar.enableSave(False)
        # statusbar = self.CreateStatusBar(2, 0)

        # NOTEBOOK
        self.notebook = wx.Notebook(self)
        self.startTab = StartTab(self.notebook)
        self.fragmentTab = FragmentTab(self.notebook, index=1)
        self.predictTab = PredictTab(self.notebook)
        self.notebook.AddPage(self.startTab, 'START Panel')
        self.notebook.AddPage(self.fragmentTab, self.fragmentTab.text)
        self.notebook.AddPage(self.predictTab, self.predictTab.text)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.EXPAND)
        sizer.Layout()
        self.SetSizer(sizer)

        self.Centre()
        # self.Show()

        # MENU EVENTS
        self.Bind(wx.EVT_MENU_OPEN, self.OnMenu)
        self.Bind(wx.EVT_MENU, self.OnLoadSARTrain, id=1)  # File
        self.Bind(wx.EVT_MENU, self.OnLoadQSARTrain, id=0)  # File
        self.Bind(wx.EVT_MENU, self.OnLoadTest, id=2)
        self.Bind(wx.EVT_MENU, self.OnLoadModel, id=wx.ID_OPEN)
        self.Bind(wx.EVT_MENU, self.OnQuit, id=wx.ID_EXIT)
        self.Bind(wx.EVT_MENU, self.OnSaveModel, id=wx.ID_SAVE)  # Save
        self.Bind(wx.EVT_MENU, self.OnSavePredEvent, id=5)
        self.Bind(wx.EVT_MENU, self.OnSavePredEvent, id=6)
        self.Bind(wx.EVT_MENU, self.OnSaveCSVModel, id=7)
        self.Bind(wx.EVT_MENU, self.OnSaveCSVCluster, id=9)
        self.Bind(wx.EVT_MENU, self.OnExportDB, id=8)
        self.Bind(wx.EVT_MENU, self.OnChangeLibrary, id=10)  # OpenBabel
        self.Bind(wx.EVT_MENU, self.OnChangeLibrary, id=11)  # RDKit
        self.Bind(wx.EVT_MENU, self.OnChangeSpeed_up, id=12)
        self.Bind(wx.EVT_MENU, OnSARpyManual, id=22)
        self.Bind(wx.EVT_MENU, OnQSARpyManual, id=20)
        self.Bind(wx.EVT_MENU, OnAbout, id=21)
        self.Bind(wx.EVT_MENU, OnDQSAR, id=30)
        # Start Panel
        self.Bind(wx.EVT_BUTTON, self.OnLoadSARTrain, id=11)
        self.Bind(wx.EVT_BUTTON, self.OnLoadQSARTrain, id=12)
        self.Bind(wx.EVT_BUTTON, self.OnStart, id=13)
        self.Bind(wx.EVT_BUTTON, self.OnStart, id=14)
        # SHORTCUT BUTTONS
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged, self.notebook)
        # self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.onPageChanging)
        self.Bind(wx.EVT_BUTTON, self.OnRefreshMenuBar, id=wx.ID_REFRESH)
        self.Bind(wx.adv.EVT_HYPERLINK, self.OnRefreshMenuBar, id=wx.ID_REFRESH)
        self.Bind(wx.EVT_BUTTON, self.OnSaveLog, id=99)  # Console.button_saveLog  > put it there!

    def OnChangeLibrary(self, event):
        if showMessage("Changing chemical library will reset all tabs. \n\n  Continue?",
                       style=wx.ICON_EXCLAMATION) == wx.ID_OK:
            self.goToPage(self.fragmentTab)  # to remove extractTab and GeneralizeTab
            self.startTab.resetTab()
            self.fragmentTab.resetTab(clear=True)
            self.predictTab.resetTab(clear=True)
            self.notebook.ChangeSelection(0)
            if event.GetId() == 10:
                chem.setLibrary('OB')
            else:
                chem.setLibrary('RDKit')

    def updateOptions(self, parameters):  # called when a saved model is opened
        if parameters['chem_library'] == 'OB':
            menu_item = self.MenuBar.FindItemById(10)
            if not menu_item.IsChecked():
                showMessage("Changing chemical library to Open Babel", style=wx.ICON_INFORMATION)
            menu_item.Check()
        else:
            menu_item = self.MenuBar.FindItemById(11)
            if not menu_item.IsChecked():
                showMessage("Changing chemical library to RDKit", style=wx.ICON_INFORMATION)
            menu_item.Check()
        self.MenuBar.FindItemById(12).Check(parameters['speed_up'])

    @staticmethod
    def OnChangeSpeed_up(event):
        showMessage("This change will take place at the next training set fragmentation.", style=wx.ICON_INFORMATION)

    def isSpeedUp(self):
        return self.MenuBar.FindItemById(12).IsChecked()

    def getModel(self):
        try:
            return self.notebook.GetCurrentPage().model
        except AttributeError:
            return None

    # def onPageChanging(self, event):
    #     event.Veto()  # disable tab clicking

    def OnPageChanged(self, event):
        self.OnRefreshMenuBar(event)
        try:
            self.notebook.GetCurrentPage().button.SetFocus()
        except:
            pass

    def OnRefreshMenuBar(self, event):
        page = self.notebook.GetCurrentPage()
        model = self.getModel()
        if model and model.ruleset:
            if page == self.predictTab:
                self.MenuBar.display('test', model.type)
            else:
                self.MenuBar.display('train', model.type)
            self.MenuBar.enableSave()
        else:
            self.MenuBar.enableSave(False)
        if page == self.predictTab:  # DISABLE 'Back' button if model changed
            prev_page = self.notebook.GetPage(self.notebook.GetSelection() - 1)
            if prev_page.model != page.model:
                page.button_left.Disable()
        self.MenuBar.Refresh()

    def OnMenu(self, event):  # 'Save menu ONLY!
        menu = event.GetMenu()
        model = self.getModel()
        if menu == self.MenuBar.delta_menu:  # D-QSAR menu
            if len(list(Path(database_path).rglob('*'))) != 0:  # it is NOT empty
                self.GetMenuBar().Enable(id=30, enable=True)  # ENABLE D-QSAR launcher
            else:
                self.GetMenuBar().Enable(id=30, enable=False)
        elif menu == self.MenuBar.TR_save_menu or menu == self.MenuBar.SAR_TR_save_menu:  # Output menu (Extractor)
            try:
                # model.trainingset.testset
                model.trainingset.predictions
            except:
                self.GetMenuBar().Enable(id=5, enable=False)
            else:
                self.GetMenuBar().Enable(id=5, enable=True)  # ENABLE save train predictions
        elif menu == self.MenuBar.TE_save_menu or menu == self.MenuBar.SAR_TE_save_menu:  # Output menu (Predictor)
            if self.predictTab.testset is not None and self.predictTab.testset.predictions:
                self.GetMenuBar().Enable(id=6, enable=True)  # ENABLE Save Predictions (Test)
            else:
                self.GetMenuBar().Enable(id=6, enable=False)
        else:
            return
        if model and model.type == 'QSAR':  # ENABLE/DISABLE save clustered model sheet
            try:
                model.tolerance
            except:
                self.GetMenuBar().Enable(id=9, enable=False)
            else:
                self.GetMenuBar().Enable(id=9, enable=True)

    def OnLoadSARTrain(self, event):
        if self.fragmentTab.loadPath():
            self.startTab.enableSAR()
            if event.GetId() == 1:  # from menu (NOT start page)
                self.OnStart(event)

    def OnLoadQSARTrain(self, event):
        if self.fragmentTab.loadPath():
            self.startTab.enableQSAR()
            if event.GetId() == 0:  # from menu (NOT start page)
                self.OnStart(event)

    def OnLoadTest(self, event):
        if self.predictTab.loadPath():
            self.startTab.enableQSAR()
            if event.GetId() == 2:
                self.goToPage(self.predictTab)
            try:
                self.GetMenuBar().Enable(id=6, enable=False)  # disable Save predictions
            except:
                pass

    def OnLoadModel(self, event):
        self.goToPage(self.predictTab)
        self.predictTab.OnLink(event)
        self.OnRefreshMenuBar(event)

    # NAVIGATION
    def goToPage(self, panel):  # ONLY for fragmentTab and predictTab
        assert panel == self.fragmentTab or panel == self.predictTab
        index = self.notebook.FindPage(panel)
        if index > 0:
            self.notebook.ChangeSelection(index)
        else:  # not found
            index = panel.index
            self.notebook.RemovePage(index)
            self.notebook.InsertPage(index, panel, panel.text)
            panel.next_page = None
            self.notebook.ChangeSelection(index)

    def OnStart(self, event):
        if event.GetId() in (0, 14):
            self.fragmentTab.setType('QSAR')
        elif event.GetId() in (1, 13):
            self.fragmentTab.setType('SAR')
        self.goToPage(self.fragmentTab)
        self.fragmentTab.widget.comboboxes[0].SetFocus()
        self.MenuBar.enableSave(False)

    def OnSaveModel(self, event=None):
        with SaveDialog(self, 'ruleset') as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                wx.BeginBusyCursor()
                path = dlg.GetPath()
                try:
                    shell.saveModel(self.getModel(), path)  # Save CURRENT page Model
                except Exception as e:
                    showMessage(dlg.msgErrorExport, exception=e, style=wx.ICON_ERROR)
                    raise
                else:
                    showMessage(dlg.msgOkExport % path.name, caption='Done!')
                finally:
                    wx.EndBusyCursor()

    def OnSaveCSVModel(self, event):
        with SaveDialog(self, 'CSV') as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                try:
                    output.output_model(self.getModel(), path)
                except Exception as e:
                    showMessage(dlg.msgErrorExport, exception=e, style=wx.ICON_ERROR)
                    raise
                else:
                    showMessage(dlg.msgOkExport % path.name, caption='Done!')

    def OnSaveCSVCluster(self, event):
        with SaveDialog(self, 'CSV') as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                try:
                    output.output_cluster_model(self.getModel(), path)
                except Exception as e:
                    showMessage(dlg.msgErrorExport, exception=e, style=wx.ICON_ERROR)
                    raise
                else:
                    showMessage(dlg.msgOkExport % path.name, caption='Done!')

    # def OnSaveCSVAnalysis(self, event):
    #     with SaveDialog(self, 'CSV') as dlg:
    #         if dlg.ShowModal() == wx.ID_OK:
    #             path = dlg.GetPath()
    #             try:
    #                 output.output_modAnalysis(self.getModel(), path)
    #             except TypeError:
    #                 showMessage("'ID' header MUST be set to save the Training Analysis.", style=wx.ICON_EXCLAMATION)
    #             except Exception as e:
    #                 showMessage(dlg.msgErrorExport, exception=e, style=wx.ICON_ERROR)
    #                 raise
    #             else:
    #                 showMessage(dlg.msgOkExport % path.name, caption='Done!')

    def OnSavePredEvent(self, event):
        eventId = event.GetId()
        if eventId == 5:
            self.savePredictions(self.notebook.GetCurrentPage().model)
        else:
            self.savePredictions(self.notebook.GetCurrentPage().model, self.predictTab.testset)

    def savePredictions(self, model, testset=None):
        dataset = testset or model.trainingset
        with SaveDialog(self, 'CSV') as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                try:
                    output.output_predictions(dataset, path, type=model.type)
                except Exception as e:
                    showMessage(dlg.msgErrorExport, exception=e, style=wx.ICON_ERROR)
                    raise
                else:
                    showMessage(dlg.msgOkExport % path.name, caption='Done!')

    def OnExportDB(self, event=None, model=None):
        model = model or self.getModel()
        if model.deep_learning:
            showMessage("Exporting a database is disabled if 'Deep learning' is enabled",
                        caption=f"{delta}-QSAR: Invalid operation", style=wx.ICON_EXCLAMATION)
            return
        if getattr(model, 'tolerance', None) is not None:
            showMessage("Exporting a database is disabled if 'Clustered View' is enabled",
                        caption=f"{delta}-QSAR: Invalid operation", style=wx.ICON_EXCLAMATION)
            return
        wait = wx.BusyCursor()
        with DBDialog(self) as param_dlg:
            if param_dlg.ShowModal() == wx.ID_OK:
                with SaveDialog(self, 'database') as dlg:
                    if dlg.ShowModal() == wx.ID_OK:
                        path = dlg.GetPath()
                        path.unlink(missing_ok=True)  # Otherwise, it ADDS rows if already exist
                        try:
                            shell.exportDB(path, model, param_dlg.getdata())
                        # except TypeError:
                        #     showMessage("'ID' header MUST be set to save the Training Analysis.",
                        #     caption=f"{delta}-QSAR: Invalid operation", style=wx.ICON_EXCLAMATION)
                        except Exception as e:
                            print("\n (Database export aborted)")
                            try:
                                path.unlink()  # does it work in windows?
                            except:
                                pass
                            showMessage(dlg.msgErrorExport, exception=e, style=wx.ICON_ERROR)
                        else:
                            showMessage(dlg.msgOkExport % path.name, caption='Done!')
        del wait

    def OnSaveLog(self, event):  # from Console
        with SaveDialog(self, 'text') as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                try:
                    with open(path, 'w', encoding='utf8') as f:
                        f.write(event.GetEventObject().GetParent().getLog())
                except Exception as e:
                    showMessage(dlg.msgErrorExport, exception=e, style=wx.ICON_ERROR)
                    raise
                else:
                    showMessage(dlg.msgOkExport % path.name, caption='Done!')

    def OnQuit(self, event):
        self.Close()


def OnSARpyManual(event):
    showMessage("Sorry, there's no SARpy manual yet!", style=wx.ICON_INFORMATION)


def OnQSARpyManual(event):
    filepath = Path(manual_html_path)
    try:
        if not filepath.exists():
            raise
        #try:
        #  if platform.system() == 'Windows':  # Windows
        #    subprocess.call(('start', filepath), shell=True)
        #except:
        import webbrowser
        webbrowser.open(str(filepath.resolve()))
    except:
        filepath = Path(manual_pdf_path)
        try:
            if not filepath.exists():
                raise
            if platform.system() == 'Darwin':  # macOS
                subprocess.call(('open', filepath))
            elif platform.system() == 'Windows':  # Windows
                subprocess.call(('start', filepath), shell=True)
                # os.startfile(filepath)
            else:  # linux variants
                subprocess.call(('xdg-open', filepath))
        except:
            showMessage("For user's manual see QSARpy_Manual.pdf", style=wx.ICON_INFORMATION)


def OnAbout(event):
    about = wx.adv.AboutDialogInfo()
    about.SetIcon(wx.Icon('./pics/SARpy_logo_256.png', wx.BITMAP_TYPE_PNG))
    about.SetName(f'{delta}-QSAR')
    about.SetVersion('0.2')
    about.SetDescription('A knowledge discovery suite for automatic SAR/QSAR rules induction.')
    about.SetCopyright('SARpy (C) 2011 - 2023 Thomas Ferrari\n'
                       'QSARpy (C) 2014 - 2023 Thomas Ferrari\n'
                       f'{delta}-QSAR (C) 2021 - 2023 Thomas Ferrari')
    # about.SetWebSite('')
    # about.SetLicense(open('./COPYING.txt', 'r').read())
    about.SetLicense("This program comes with ABSOLUTELY NO WARRANTY.\n"
                     "This is free software, and you are welcome to redistribute it under certain conditions.\n"
                     "See 'COPYING.txt' for details.\n")
    about.AddDeveloper('Thomas Ferrari\n Email: thedataconspiracy@gmail.com\n')
    # about.AddDeveloper('thedataconspiracy@gmail.com')
    about.AddArtist('Vittorio Castiglioni')
    about.AddArtist('Dalila Facchinetti\n')
    wx.adv.AboutBox(about)


def launchDeltaQSAR():
    try:
        subprocess.call(Path(DQSAR_path).resolve())
    except Exception as e:
        showMessage("Error: try to execute the file manually:\n\n" + DQSAR_path,
                    exception=e, style=wx.ICON_ERROR)


def OnDQSAR(event):
    SubThread(target=launchDeltaQSAR).start()


def close():
    return wx.ID_CANCEL


def showMessage(text, caption=None, exception=None, style=wx.OK):
    if caption is None:
        caption = ""
        if style == wx.ICON_INFORMATION:
            caption = f"{delta}-QSAR: Information"
        elif style == wx.ICON_EXCLAMATION:
            caption = f"{delta}-QSAR: Warning"
        elif style == wx.ICON_ERROR:
            caption = f"{delta}-QSAR: Error!"
    if exception:
        text += '\n\n Exception: ' + repr(exception)
    msg_dlg = wx.MessageDialog(None, text, caption, style)
    msg_dlg.Bind(wx.EVT_CLOSE, close)
    return msg_dlg.ShowModal()
    # return wx.MessageDialog(None, text, caption, style).ShowModal()


# PROGRESS BAR NEEDS DIFFERENT THREADS
# progress = wx.ProgressDialog("Exporting database", "Please wait...", maximum=len(model.neuron_log), parent=self)
#      parent=self, style=wx.PD_SMOOTH | wx.PD_AUTO_HIDE | wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
# progress.ShowModal()
# for i in delta_DB.exportDB(model, path, parameters):
#     print(i)
#     progress.Update(i)
#     wx.Yield()
# progress.Destroy()

# def start(func, *args):  # helper method to run a function in another thread
#     thread = threading.Thread(target=func, args=args)
#     thread.setDaemon(True)
#     thread.start()



################################################################################


# MySplash = wx.adv.SplashScreen(bmp_splash, wx.adv.SPLASH_CENTRE_ON_SCREEN | wx.adv.SPLASH_TIMEOUT, 3000, mainFrame)
# wx.SafeYield()


class Splash(SPLASH.AdvancedSplash):

    def __init__(self, parent):
        img_DQSAR_splash.ConvertAlphaToMask()
        bmp_splash = wx.Bitmap(img_DQSAR_splash)
        super().__init__(parent, bitmap=bmp_splash, timeout=3000)

    def OnPaint(self, event):
        # self.GetParent().Hide()
        super().OnPaint(event)

    def OnNotify(self, event):
        super().OnNotify(event)
        self.GetParent().Show()


def unblockDlls():
    dlls = glob(r"./*.dll")
    for dll in dlls:
        try:
            os.remove(dll + ':Zone.Identifier')
        except FileNotFoundError:  # The ADS did not exist or already unblocked
            pass


if __name__ == '__main__':
    unblockDlls()
    chem.setLibrary()  # OpenBabel or RDKit (default='OB', otherwise RDKit)
    app = wx.App()
    fonts = MyFonts()

    mainFrame = QSARpyFrame(None, title=f'{delta}-QSAR')
    mainFrame.SetMinSize(frame_size)
    splash = Splash(mainFrame)
    app.MainLoop()
