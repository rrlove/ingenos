from collections import namedtuple

Inversion = namedtuple('Inversion',['arm','proximal_start','proximal_end','distal_start','distal_end'])
inversionDict = {}
offset = 500

inversionDict["2La"] = Inversion(b'2L',20524058,20528089,42165182,42165532)
inversionDict["2Rb"] = Inversion(b'2R',19023925,19027916,26747166,26758676)
inversionDict["2Rc"] = Inversion(b'2R',26750000,26784943,31450000,31473100)
inversionDict["2Ru"] = Inversion(b'2R',31473000,31483751,35504441,35505236)
inversionDict["2Rd"] = Inversion(b'2R',31480000-offset,31480000+offset,42600000-offset,42600000+offset)
inversionDict["2Rj"] = Inversion(b'2R',3262186,3262296,15750716,15750717)
