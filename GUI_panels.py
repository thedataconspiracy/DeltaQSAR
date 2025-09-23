#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021  Thomas Ferrari


import wx
import platform
import shell
from math import log, exp

# import warnings
# warnings.filterwarnings("ignore")

DECIMAL = shell.DECIMAL

combobox_minsize = (150, -1)
spinctrl_minsize = (80, -1)
button_minsize = (-1, 50)
indentation = (5, -1)
wrap = 500
app = wx.App()


class MyFonts:

    def __init__(self):
        ps = 11
        self.normal = wx.Font(wx.FontInfo(ps))
        self.bold = wx.Font(wx.FontInfo(ps).Bold())
        self.italic = wx.Font(wx.FontInfo(ps).Slant())
        self.smaller_italic = wx.Font(wx.FontInfo(ps - 1).Slant())
        self.smaller = wx.Font(wx.FontInfo(ps - 1))
        self.bigger = wx.Font(wx.FontInfo(ps + 1))
        self.title = wx.Font(wx.FontInfo(ps + 3).Bold())
        self.boldino = wx.Font(wx.FontInfo(ps - 1).Bold())


fonts = MyFonts()


class StartPanel(wx.Panel):

    SAR = {'intro': "Submit a training set to extract Structural Alerts for the endpoint of interest.",
           'footer_1': "ordered list of SMARTS",
           'footer_2': "with the relative endpoint class."}
    QSAR = {'intro': "Submit a training set to extract a model for the endpoint of interest as a set of structural rules.",
            'footer_1': "set of Base Structures (with activity values)",
            'footer_2': "set of Modulators (with \u03B4-activity values)"}

    def __init__(self, parent, title):
        super().__init__(parent)
        self.box = wx.StaticBox(self, -1, title)
        self.st_inputstitle = wx.StaticText(self, -1, "INPUT")
        self.st_outputstitle = wx.StaticText(self, -1, "OUTPUT")
        # self.st_ruleset = wx.StaticText(self, -1, "Ruleset")
        self.st_arrow_in = wx.StaticText(self, -1, "\u2193")
        self.st_arrow_out = wx.StaticText(self, -1, "\u2193")
        bmp = wx.Bitmap('./pics/dataset.png', wx.BITMAP_TYPE_ANY)

        if 'QSAR' in title:
            self.button_load = wx.BitmapButton(self, 12, bmp, size=(bmp.GetWidth() + 10, bmp.GetHeight() + 10))
            # self.mode = 'test'
            self.st_intro = wx.StaticText(self, -1, self.QSAR['intro'])
            self.st_input = wx.StaticText(self, -1, "Dataset")
            self.st_output = wx.StaticText(self, -1, "Ruleset")
            self.button = wx.Button(self, 14, "QSAR")
            self.st_footertitle = wx.StaticText(self, -1, 'Ruleset :')
            self.st_footer1 = wx.StaticText(self, -1, self.QSAR['footer_1'])
            self.st_footer2 = wx.StaticText(self, -1, self.QSAR['footer_2'])
            # self.st_footertitle = wx.StaticText(self, -1, 'Prediction range :')
            # self.st_footer1 = wx.StaticText(self, -1, 'Any structure S such that:')
            # self.st_footer2 = wx.StaticText(self, -1, 'S = Base Structure \u00B1 Modulator(s)')
        elif 'SAR' in title:
            self.button_load = wx.BitmapButton(self, 11, bmp, size=(bmp.GetWidth() + 10, bmp.GetHeight() + 10))
            # self.mode = 'train'
            self.st_intro = wx.StaticText(self, -1, self.SAR['intro'])
            # self.st_input = wx.StaticText(self, -1, "Dataset            ")
            self.st_input = wx.StaticText(self, -1, "Dataset")
            self.st_output = wx.StaticText(self, -1, "Ruleset")
            self.button = wx.Button(self, 13, "SAR")
            self.st_footertitle = wx.StaticText(self, -1, 'Ruleset :')
            self.st_footer1 = wx.StaticText(self, -1, self.SAR['footer_1'])
            self.st_footer2 = wx.StaticText(self, -1, self.SAR['footer_2'])

        self.button_load.SetToolTip(wx.ToolTip("LOAD dataset"))
        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.st_intro.SetFont(fonts.normal)
        self.st_inputstitle.SetFont(fonts.boldino)
        self.st_outputstitle.SetFont(fonts.boldino)
        self.st_input.SetFont(fonts.normal)
        self.st_output.SetFont(fonts.normal)
        self.st_arrow_in.SetFont(fonts.bigger)
        self.st_arrow_out.SetFont(fonts.bigger)
        self.button.SetMinSize(button_minsize)
        self.button.Disable()
        self.st_footertitle.SetFont(fonts.boldino)
        self.box.SetFont(fonts.bold)
        self.st_footer2.SetFont(fonts.smaller_italic)
        self.st_footer1.SetFont(fonts.smaller_italic)
        self.st_intro.Wrap(wrap)

    def __do_layout(self):
        sbs = wx.StaticBoxSizer(self.box, wx.VERTICAL)
        sbs.Add(self.st_intro, 1, wx.EXPAND | wx.ALL, 15)

        sbs.Add((-1, 1), 1)  # expandable spacer

        vbox_in = wx.BoxSizer(wx.VERTICAL)
        vbox_in.Add(self.st_inputstitle, 0, wx.ALIGN_LEFT)
        vbox_in.Add(self.st_input, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        vbox_load = wx.BoxSizer(wx.VERTICAL)
        vbox_load.Add(self.button_load, 0, wx.ALIGN_LEFT)
        hbox_input = wx.BoxSizer(wx.HORIZONTAL)
        hbox_input.Add((1, -1), 1)
        hbox_input.Add(vbox_load, 1)
        hbox_input.Add(vbox_in, 0, wx.ALIGN_BOTTOM | wx.BOTTOM, 10)
        hbox_input.Add((1, -1), 1)
        sbs.Add(hbox_input, 0, wx.EXPAND)
        sbs.Add(self.st_arrow_in, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 10)

        hbox_buttons = wx.BoxSizer(wx.HORIZONTAL)
        hbox_buttons.Add((1, -1), 1)
        hbox_buttons.Add(self.button, 2, flag=wx.EXPAND)
        hbox_buttons.Add((1, -1), 1)
        sbs.Add(hbox_buttons, 2, flag=wx.EXPAND)

        sbs.Add(self.st_arrow_out, 0, wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 10)
        hbox_output = wx.BoxSizer(wx.HORIZONTAL)
        vbox_out = wx.BoxSizer(wx.VERTICAL)
        vbox_out.Add(self.st_outputstitle, 0, wx.ALIGN_RIGHT)
        vbox_out.Add(self.st_output, 0, wx.ALIGN_RIGHT | wx.TOP, 5)
        hbox_output.Add(vbox_out, 1)
        hbox_output.Add((1, -1), 1)
        sbs.Add(hbox_output, 0, wx.EXPAND)

        sbs.Add((-1, 1), 1)  # expandable spacer

        vbox_footer = wx.BoxSizer(wx.VERTICAL)
        vbox_footer.Add(self.st_footer1)
        vbox_footer.Add(self.st_footer2)
        hbox_footer = wx.BoxSizer(wx.HORIZONTAL)
        hbox_footer.Add(self.st_footertitle)
        hbox_footer.Add(vbox_footer, 0, wx.LEFT, 10)
        sbs.Add(hbox_footer, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 20)
        sbs.Layout()
        self.SetSizer(sbs)

    def enable(self, enable=True):
        self.button.Enable(enable)
        self.button.SetFocus()


class Console(wx.Panel):
    labels = {'path_label': "PATH :",
              'console_text': "\n Load a Data set from the File menu to start."}
    blank = True

    def __init__(self, parent, labels={}):
        super().__init__(parent)
        self.path = None
        self.font = None
        self.labels.update(labels)
        self.st_dataset = wx.StaticText(self, -1, self.labels['path_label'])
        self.tc_path = wx.TextCtrl(self, -1, '', style=wx.TE_READONLY)
        self.st_consoletitle = wx.StaticText(self, -1, "CONSOLE OUTPUT")
        self.tc_logArea = wx.TextCtrl(self, -1, self.labels['console_text'],
                                      style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH)

        bmp_clear = wx.Bitmap('./pics/clear.png', wx.BITMAP_TYPE_ANY)
        bmp_log = wx.Bitmap('./pics/log1.png', wx.BITMAP_TYPE_ANY)
        self.button_clear = wx.BitmapButton(self, -1, bmp_clear)
        self.button_saveLog = wx.BitmapButton(self, 99, bmp_log)
        self.button_clear.SetToolTip(wx.ToolTip("CLEAR console"))
        self.button_saveLog.SetToolTip(wx.ToolTip("SAVE log file"))
        # self.button_clear = wx.Button(self, -1, "Clear")
        # self.button_saveLog = wx.Button(self, 99, "Save LOG")

        self.button_clear.Bind(wx.EVT_BUTTON, self.OnClear)

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.tc_path.Disable()
        self.st_consoletitle.SetFont(fonts.bold)
        self.st_dataset.SetFont(fonts.smaller)
        # In Win Courier New is ugly, in Linux Lucida Sans Typewriter is not monospaced
        if platform.system() == 'Windows':
            facename = 'Lucida Sans Typewriter'
        else:
            facename = 'Courier New'
        self.font = wx.Font(wx.FontInfo(10).Family(wx.FONTFAMILY_TELETYPE).FaceName(facename))
        self.tc_logArea.SetFont(self.font)
        # self.tc_logArea.SetBackgroundColour(wx.Colour(0, 0, 0))
        # self.tc_logArea.SetForegroundColour(wx.Colour(255, 255, 255))
        self.button_clear.SetMinSize(button_minsize)
        self.button_saveLog.SetMinSize(button_minsize)
        self.button_clear.Disable()
        self.button_saveLog.Disable()

    def __do_layout(self):
        vbox_dx = wx.BoxSizer(wx.VERTICAL)
        hbox_inputs = wx.BoxSizer(wx.HORIZONTAL)
        hbox_inputs.Add(self.st_dataset, 0, wx.ALIGN_CENTER_VERTICAL)
        hbox_inputs.Add(self.tc_path, 1, wx.LEFT, 10)
        vbox_dx.Add(hbox_inputs, 0, wx.EXPAND)
        vbox_dx.Add((-1, 15))

        vbox_dx.Add(self.st_consoletitle, 0)
        vbox_dx.Add(self.tc_logArea, 1, wx.EXPAND | wx.TOP, border=5)
        vbox_dx.Add((-1, 10))

        hbox_buttons = wx.BoxSizer(wx.HORIZONTAL)
        hbox_buttons.Add(self.button_clear, 1, wx.EXPAND)
        hbox_buttons.Add(self.button_saveLog, 1, wx.EXPAND)
        vbox_dx.Add(hbox_buttons, 0, wx.EXPAND)
        vbox_dx.Layout()
        self.SetSizer(vbox_dx)

    def setStatus(self, text):
        self.tc_path.SetValue(text)
        # self.tc_path.SetInsertionPointEnd()
        self.tc_path.Enable()

    def write(self, text, color=None, bold=False):
        if self.blank:
            self.tc_logArea.Clear()
            self.blank = False
        self.tc_logArea.SetDefaultStyle(wx.TextAttr(wx.BLACK, wx.NullColour, self.font))
        if color:
            self.tc_logArea.SetDefaultStyle(wx.TextAttr(color))
        if bold:
            self.tc_logArea.SetDefaultStyle(wx.TextAttr(wx.BLACK, wx.NullColour, self.font.Bold()))
        # self.tc_logArea.write(text)
        wx.CallAfter(self.tc_logArea.write, text)
        self.tc_logArea.SetFocus()
        self.tc_logArea.SetInsertionPointEnd()
        self.tc_logArea.Refresh()
        wx.SafeYield()  # to flush
        self.button_clear.Enable()
        self.button_saveLog.Enable()

    def flush(self):
        self.tc_logArea.flush()  # stdout redirected here

    def getLog(self):
        return self.tc_logArea.GetValue()

    def OnClear(self, event=None, prompt=None):
        self.tc_logArea.Clear()
        self.tc_path.Clear()
        self.tc_path.Disable()
        if prompt is not None:
            self.write(prompt)
        self.blank = True
        self.button_clear.Disable()
        self.button_saveLog.Disable()


class HeaderPanel(wx.Panel):
    cb_labels = ["Select the CAS header (or any ID)...",
                 "Select the SMILES header...",
                 "Select the ENDPOINT header..."]

    def __init__(self, parent, title, outlier_checkbox=False):
        super().__init__(parent)
        self.st_headerstitle = wx.StaticText(self, -1, title)
        self.combobox_ID = wx.ComboBox(self, -1, choices=[], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.combobox_SMILES = wx.ComboBox(self, -1, choices=[], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.combobox_target = wx.ComboBox(self, -1, choices=[], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.st_starId = wx.StaticText(self, -1, '*')
        self.st_starSmiles = wx.StaticText(self, -1, '*')
        self.st_starTarget = wx.StaticText(self, -1, '*')
        self.st_requiredFields = wx.StaticText(self, -1, '* Required fields')
        # self.checkbox_outliers = wx.CheckBox(self, -1, label="Preprocess outliers")
        self.button_load = wx.Button(self, -1, "Load")
        self.comboboxes = (self.combobox_ID, self.combobox_SMILES, self.combobox_target)

        # PROPERTIES
        self.st_headerstitle.SetFont(fonts.bold)
        self.st_requiredFields.SetFont(fonts.smaller)
        # self.checkbox_outliers.SetValue(False)
        # self.checkbox_outliers.SetMinSize((-1, 50))

        # EVENTS
        for cb in self.comboboxes:
            cb.Bind(wx.EVT_COMBOBOX, self.OnCBevent)
        # self.checkbox_outliers.Bind(wx.EVT_CHECKBOX, self.OnCBevent)

        # if not outlier_checkbox:
        #     self.checkbox_outliers.Hide()
        self.reset()

        hbox_id = wx.BoxSizer(wx.HORIZONTAL)
        hbox_id.Add(indentation)
        hbox_id.Add(self.combobox_ID, 1)
        hbox_id.Add(self.st_starId,
                    flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.RESERVE_SPACE_EVEN_IF_HIDDEN, border=10)

        hbox_smiles = wx.BoxSizer(wx.HORIZONTAL)
        hbox_smiles.Add(indentation)
        hbox_smiles.Add(self.combobox_SMILES, 1)
        hbox_smiles.Add(self.st_starSmiles,
                        flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.RESERVE_SPACE_EVEN_IF_HIDDEN, border=10)

        hbox_activity = wx.BoxSizer(wx.HORIZONTAL)
        hbox_activity.Add(indentation)
        hbox_activity.Add(self.combobox_target, 1)
        hbox_activity.Add(self.st_starTarget,
                          flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.RESERVE_SPACE_EVEN_IF_HIDDEN, border=10)

        vbox_load = wx.BoxSizer(wx.VERTICAL)
        # vbox_load.Add(self.st_requiredFields, flag=wx.ALIGN_RIGHT)
        vbox_load.Add((-1, 10))
        vbox_load.Add(self.button_load, flag=wx.ALIGN_RIGHT)
        hbox_button = wx.BoxSizer(wx.HORIZONTAL)
        hbox_button.Add((10, -1))
        # hbox_button.Add(self.checkbox_outliers)
        hbox_button.Add(self.st_requiredFields)
        hbox_button.Add(vbox_load, 1)
        hbox_button.Add((20, -1))

        vbox_headers = wx.BoxSizer(wx.VERTICAL)
        vbox_headers.Add((-1, 5))
        vbox_headers.Add(self.st_headerstitle, flag=wx.TOP, border=10)
        vbox_headers.Add((-1, 5))
        vbox_headers.Add(hbox_id, flag=wx.TOP | wx.EXPAND, border=10)
        vbox_headers.Add(hbox_smiles, flag=wx.TOP | wx.EXPAND, border=10)
        vbox_headers.Add(hbox_activity, flag=wx.TOP | wx.EXPAND, border=10)
        vbox_headers.Add((-1, 10))
        vbox_headers.Add(hbox_button, flag=wx.EXPAND)
        vbox_headers.Layout()
        self.SetSizer(vbox_headers)

    def setRequirements(self, ID=True, smiles=True, target=True):
        self.st_starId.Show(ID)
        self.combobox_ID.isRequired = ID
        self.st_starSmiles.Show(smiles)
        self.combobox_SMILES.isRequired = smiles
        self.st_starTarget.Show(target)
        self.combobox_target.isRequired = target
        self.combobox_ID.SetFocus()
        if not smiles:
            self.combobox_SMILES.Disable()
        if not (ID or smiles or target):
            self._enableButton()

    def OnCBevent(self, event):
        for cb in self.comboboxes:
            if cb.GetSelection() == 0 and cb.isRequired:
                self.button_load.Disable()
                cb.SetFocus()
                return
        self._enableButton()

    def _enableButton(self):
        self.button_load.Enable()
        self.button_load.SetFocus()

    def setHeaders(self, headers):
        for cb in self.comboboxes:
            cb.AppendItems(headers)
        self.enable()

    def getHeaders(self, index=None):
        idKey = smilesKey = targetKey = None
        if index:
            if self.comboboxes[index].GetSelection() != 0:
                return self.comboboxes[index].GetStringSelection()
            else:
                return None
        else:
            if self.combobox_ID.GetSelection() != 0:
                idKey = self.combobox_ID.GetStringSelection()
            if self.combobox_SMILES.GetSelection() != 0:
                smilesKey = self.combobox_SMILES.GetStringSelection()
            if self.combobox_target.GetSelection() != 0:
                targetKey = self.combobox_target.GetStringSelection()
        return idKey, smilesKey, targetKey

    def reset(self):
        for i, cb in enumerate(self.comboboxes):
            cb.Clear()
            cb.Append(self.cb_labels[i])
            cb.SetSelection(0)
        self.st_starId.Hide()
        self.st_starSmiles.Hide()
        self.st_starTarget.Hide()
        # self.checkbox_outliers.SetValue(False)
        self.disable()

    def disable(self):
        for cb in self.comboboxes:
            cb.Disable()
        self.button_load.Disable()
        # self.checkbox_outliers.Disable()

    def enable(self):
        for cb in self.comboboxes:
            cb.Enable()
        # self.checkbox_outliers.Enable()


class SAR_ExtractPanel(wx.Panel):

    labels1 = {'left_label': "MAX ERROR :",
               'right_label': "",
               'left_label_2': "(precision)",
               'right_label_2': "(coverage)",
               'value': "%"
               }
    labels1_mono = {'left_label': "MAX ERROR :",
               'right_label': "",
               'left_label_2': "(specificity)",
               'right_label_2': "(sensitivity)",
               'value': "%"
                   }
    labels2 = {'left_label': "precision",
               'right_label': "recall",
               'left_label_2': "",
               'right_label_2': "",
               'value': ""
               }

    def __init__(self, parent, targets):
        super().__init__(parent)
        self.cb_target = wx.ComboBox(self, -1, choices=['ALL classes'], style=wx.CB_DROPDOWN | wx.CB_READONLY)

        self.st_slider_title = wx.StaticText(self, -1, "EXTRACTION:")
        self.sliderPanel1 = SliderPanel(self, custom_labels=self.labels1)
        self.sliderPanel2 = SliderPanel(self, custom_labels=self.labels2)
        self.sliderPanel1.spin.SetIncrement(1)
        self.sliderPanel2.spin.SetIncrement(0.01)

        self.cb_target.AppendItems([f"'{target}'" for target in targets])  # adding quotes
        self.cb_target.SetSelection(0)
        self.cb_target.Bind(wx.EVT_COMBOBOX, self.OnSelectTarget)

        self.cb_target.SetMinSize(combobox_minsize)
        self.st_slider_title.SetFont(fonts.bold)
        self.reset()
        self._do_layout()

    def _do_layout(self):
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.st_slider_title, 0)
        hbox.Add((1, -1), 1)  # expandable spacer
        hbox.Add(self.cb_target, 0)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(hbox, flag=wx.EXPAND | wx.TOP, border=15)
        vbox.Add((-1, 10))
        vbox.Add(self.sliderPanel1, flag=wx.EXPAND)
        vbox.Add(self.sliderPanel2, flag=wx.EXPAND)
        vbox.Layout()
        self.SetSizer(vbox)

    def reset(self, target=None):
        if target is None:
            self.sliderPanel1.init(0, 25, 50)
        else:
            self.sliderPanel1.init(0, 50, 100)
        self.sliderPanel2.init(0, 0.05, 0.1)

    def OnSelectTarget(self, event=None):
        if self.cb_target.GetSelection() == 0:
            self.sliderPanel1.setSliderLabels(self.labels1)
            self.reset()
        else:
            self.sliderPanel1.setSliderLabels(self.labels1_mono)
            target = self.cb_target.GetStringSelection()
            self.reset(target)
        self.sliderPanel1.GetSizer().Layout()  # redo layout after changing labels

    def getDict(self):
        if self.cb_target.GetSelection() == 0:
            target = None
        else:
            target = self.cb_target.GetStringSelection()[1:-1]  # removing quotes
        max_error = round(self.sliderPanel1.getValue() / 100, DECIMAL)
        return {'mono_target': target, 'max_error': max_error, 'beta': self.sliderPanel2.getValue()}
        # return target, self.sliderPanel1.getValue(), self.sliderPanel2.getValue()


class SliderPanel(wx.Panel):
    base_labels = {'title': "",
                   'left_label': "",
                   'right_label': "",
                   'left_label_2': "",
                   'right_label_2': "",
                   'value': "Value :"
                   }
                   # 'text': ""}

    def __init__(self, parent, custom_labels={}, custom_sizer=False):
        # If custom_sizer is set, you must call _do_layout() manually
        super().__init__(parent)
        labels = dict(SliderPanel.base_labels)
        labels.update(custom_labels)
        self.st_title = wx.StaticText(self, -1, labels['title'])
        # self.st_intro = wx.StaticText(self, -1, "Increase the tolerance of the model to simplify the Ruleset")
        self.st_max1 = wx.StaticText(self, -1, labels['right_label'])
        self.st_max2 = wx.StaticText(self, -1, labels['right_label_2'])
        self.st_min1 = wx.StaticText(self, -1, labels['left_label'])
        self.st_min2 = wx.StaticText(self, -1, labels['left_label_2'])
        self.spin = wx.SpinCtrlDouble(self, inc=0.1)
        self.st_value = wx.StaticText(self, -1, labels['value'])
        # self.infinite = wx.StaticText(self, -1, "")
        self.slider = wx.Slider(self)
        # tip = wx.ToolTip("Extract only modulators with\ntraining error below the threshold.")
        # self.slider.SetToolTip(tip)

        # PROPERTIES
        self.st_title.SetFont(fonts.bold)
        # self.st_value.SetFont(fonts.smaller)
        self.st_max1.SetFont(fonts.smaller)
        self.st_max2.SetFont(fonts.smaller)
        self.st_min1.SetFont(fonts.smaller)
        self.st_min2.SetFont(fonts.smaller)
        # self.st_text.SetFont(fonts.smaller)
        self.spin.SetMinSize(spinctrl_minsize)
        self.spin.SetDigits(DECIMAL)  # set number of decimal digits
        # EVENTS
        self.slider.Bind(wx.EVT_SLIDER, self.OnSliding)
        self.spin.Bind(wx.EVT_SPINCTRLDOUBLE, self.OnSpinning)

        # if isinstance(self, OptimalSliderPanel):
        #     self.custom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # else:
        self._do_layout()

    def _do_layout(self):
        vbox_min = wx.BoxSizer(wx.VERTICAL)
        # vbox_min.Add(self.st_min1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        # vbox_min.Add(self.st_min2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        vbox_min.Add(self.st_min1, 0, wx.ALIGN_RIGHT)
        vbox_min.Add(self.st_min2, 0, wx.ALIGN_RIGHT)

        hbox_value = wx.BoxSizer(wx.HORIZONTAL)
        hbox_value.Add(self.st_value, 0, wx.ALIGN_CENTER_VERTICAL)
        # hbox_value.Add(self.infinite, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL, 20)  # infinite
        hbox_value.Add(self.spin, 0, wx.LEFT, 5)

        vbox_slider = wx.BoxSizer(wx.VERTICAL)
        vbox_slider.Add(self.slider, 0, wx.EXPAND)
        vbox_slider.Add((-1, 15))
        vbox_slider.Add(hbox_value, 0, wx.TOP | wx.ALIGN_RIGHT)  #

        vbox_max = wx.BoxSizer(wx.VERTICAL)
        # vbox_max.Add(self.st_max1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        # vbox_max.Add(self.st_max2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        vbox_max.Add(self.st_max1, 0, wx.ALIGN_LEFT)
        vbox_max.Add(self.st_max2, 0, wx.ALIGN_LEFT)

        hbox_slider = wx.BoxSizer(wx.HORIZONTAL)
        # hbox_slider.Add(indentation)
        hbox_slider.Add(vbox_min, 1)  #, wx.TOP, 5)
        hbox_slider.Add(vbox_slider, 2, wx.LEFT | wx.RIGHT, 10)
        hbox_slider.Add(vbox_max, 1)  #, wx.TOP, 5)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add((-1, 5))
        # vbox.Add(self.st_note, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)
        # vbox.Add((-1, 5))
        # vbox.Add(self.st_intro, flag=wx.LEFT, border=15)
        # vbox.Add((-1, 20))
        if self.st_title.GetLabel():
            vbox.Add(self.st_title, flag=wx.TOP, border=10)
            vbox.Add((-1, 1), 2)  # expandable spacer
        else:
            vbox.Add((-1, 1), 1)  # expandable spacer
        vbox.Add(hbox_slider, 2, flag=wx.RIGHT | wx.LEFT | wx.EXPAND, border=10)
        # vbox.Add((-1, 1), 1)  # expandable spacer
        # vbox.Add((-1, 40))
        try:
            vbox.Add(self.custom_sizer, flag=wx.ALIGN_CENTER_HORIZONTAL)
        except AttributeError:
            pass
        vbox.Layout()
        self.SetSizer(vbox)

    def setSliderLabels(self, label_dict):
        self.st_min1.SetLabel(label_dict['left_label'])
        self.st_max1.SetLabel(label_dict['right_label'])
        self.st_min2.SetLabel(label_dict['left_label_2'])
        self.st_max2.SetLabel(label_dict['right_label_2'])
        self.st_value.SetLabel(label_dict['value'])

    def init(self, start=0, default=0.5, end=1):
        # VALUES ARE ROUNDED!
        try:
            assert DECIMAL > 0  # try to set rounder value
            r_start = round(start, DECIMAL - 1)
            r_value = round(default, DECIMAL - 1)
            r_end = round(end, DECIMAL - 1)
            assert r_start <= r_value < r_end
        except AssertionError:
            r_start = round(start, DECIMAL)
            r_value = round(default, DECIMAL)
            r_end = round(end, DECIMAL)
        if end >= 50:
            self.spin.SetDigits(0)
        else:
            self.spin.SetDigits(DECIMAL)
        self.spin.SetMin(r_start)
        self.spin.SetMax(r_end)
        self.setValue(r_value)

    def setValue(self, value):
        self.spin.SetValue(value)
        self.updateSlider()

    def updateSlider(self):
        self.slider.SetValue(self.spin.GetValue() * 100 / self.spin.GetMax())
        # pass

    def getValue(self):
        return round(self.spin.GetValue(), DECIMAL)

    def OnSliding(self, event):
        output = self.slider.GetValue() * self.spin.GetMax() / 100
        self.spin.SetValue(round(output, 2))

    def OnSpinning(self, event):
        self.updateSlider()

    def enable(self, enable=True):
        self.slider.Enable(enable)
        self.spin.Enable(enable)
        # if not enable:
        #     self.slider_old_value = self.slider.GetValue()
        #     self.slider.SetValue(self.slider.GetMin())
        # else:
        #     try:
        #         self.slider.SetValue(self.slider_old_value)
        #     except AttributeError:
        #         pass


# class OptimalSliderPanel(SliderPanel):
#
#     def __init__(self, parent, custom_labels={}):
#         super().__init__(parent, custom_labels)
#         self.st_text = wx.StaticText(self, -1, "Binary classification:")
#         self.toggle_button = wx.ToggleButton(self, -1, "OPTIMAL tradeoff")
#         self.st_text.SetFont(fonts.smaller)
#         self.toggle_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnOptimal)
#         # pre-layout
#         # self.custom_sizer = wx.BoxSizer(wx.HORIZONTAL)
#         self.custom_sizer.Add(self.st_text, flag=wx.ALIGN_CENTER_VERTICAL)
#         self.custom_sizer.Add(self.toggle_button, border=200)
#         self._do_layout()
#
#     def OnOptimal(self, event):
#         self.enable(not self.toggle_button.GetValue())

    # def updateSlider(self):
    #     value = self.spin.GetValue()
    #     self.slider.SetValue(log(value) * 100 / log(self.spin.GetMax()))
    #     if value == self.spin.GetMax():
    #         self.set_inf()
    #
    # def OnSliding(self, event):
    #     value = self.slider.GetValue()
    #     output = exp(value * log(self.spin.GetMax()) / 100)
    #     self.spin.SetValue(round(output, 2))
    #     if value == self.slider.GetMax():
    #         self.set_inf()
    #     else:
    #         self.spin.Show()
    #         self.infinite.SetLabel("")
    #
    # def getValue(self):
    #     value = self.spin.GetValue()
    #     if value == self.spin.GetMax():
    #         return float('inf')
    #     return round(value, DECIMAL)
    #
    # def set_inf(self):
    #     self.spin.Hide()
    #     self.infinite.SetLabel("Infinite")


