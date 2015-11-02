Call SetLocale("en-us")
Call getDocument().beginUndoGroup("Set Default Units", true)
Call getDocument().setDefaultLengthUnit("Millimeters")
Call getDocument().endUndoGroup()
Call getDocument().setCurveSmoothnessAngle (1)

'COSTRUZIONE MEZZA CAVA

