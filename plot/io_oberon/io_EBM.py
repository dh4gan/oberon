# Written 15/1/14 Duncan Forgan
# Methods for reading EBM data

# EBM Log File dictionaries

logvariablekeys = ("t","minT", "maxT", "meanT", "meanQ","meanA","MeanIR", "meanS", 
                "meanhab", "meantidal", "a","ecc","inc","argper","longascend", "meananom", "obliq", "precession", "co2", "diffusion")
logvariablenames = ("Time (yr)", "minimum T (K)", "maximum T (K)", "mean T (K)", r"Net Heating ($erg\,s^{-1}\, cm^{-2}$)", r"Albedo", 
                 r"Mean Cooling ($erg\,s^{-1}\, cm^{-2}$)", r"Mean Insolation ($erg\,s^{-1}\, cm^{-2}$)", "Habitability Index", "Tidal Heating",
                 "Semimajor Axis", "Eccentricity", "Inclination","Argument of Periapsis", "Longitude of the Ascending Node",  "Mean Anomaly", "Obliquity", "Precession", r"$P_{CO_2}$", r"Diffusion Constant")
logvariablecolumns = range(20)


lognamedict = dict(zip(logvariablekeys,logvariablenames))
logcoldict = dict(zip(logvariablekeys,logvariablecolumns))


# EBM snapshot dictionaries

snapvariablekeys = ("x","lat", "T", "C", "Q","IR","Albedo", "S", 
                "tau", "ice","hab")
snapvariablenames = ("x", r"$\lambda$", "T (K)", "C (/)", r"Net Heating ($erg\,s^{-1}\, cm^{-2}$)",  
                 r"IR Cooling ($erg\,s^{-1}\, cm^{-2}$)",r"Albedo", r"Mean Insolation ($erg\,s^{-1}\, cm^{-2}$)", 
                 "Optical Depth", r"$f_{ice}$","Habitability Index")                  
snapvariablecolumns = (0,1,2,3,4,5,6,7,8,9,10)

snapnamedict = dict(zip(snapvariablekeys,snapvariablenames))
snapcoldict = dict(zip(snapvariablekeys,snapvariablecolumns))


def select_variables(variablekeys, coldict, labeldict):
    print "Which variable for the x-axis?"

    for i in range(len(variablekeys)):
        print variablekeys[i],":\t \t \t", labeldict[variablekeys[i]]

    keyword = raw_input("Enter appropriate keyword:   ")

    xkey = keyword
    ix = coldict[keyword]

    print "Which variable for the y-axis?"

    for i in range(len(variablekeys)):
        if(i!=ix): print variablekeys[i],":\t \t", labeldict[variablekeys[i]]

    keyword = raw_input("Enter appropriate keyword:   ")

    ykey = keyword
    iy = coldict[keyword]
        
    return xkey, ykey, ix,iy


def select_multiple_variables(variablekeys, coldict, labeldict):
    print "Which variable for the x-axis?"

    for i in range(len(variablekeys)):
        print variablekeys[i],":\t \t \t", labeldict[variablekeys[i]]

    keyword = raw_input("Enter appropriate keyword:   ")

    xkey = keyword
    ix = coldict[keyword]

    print "Which variables for the y-axis?"

    for i in range(len(variablekeys)):
        if(i!=ix): print variablekeys[i],":\t \t", labeldict[variablekeys[i]]

    keywords = raw_input("Enter keywords separated by spaces:   ")

    # Otherwise, parse keyword string into individual choices
    ykeys = keywords.split()
        
    
    iylist = []
    for word in ykeys:
        iylist.append(coldict[word])
    
    nyvar = len(iylist)        
    return xkey, ykeys, ix, iylist,nyvar


