CATEGORY_A=[]
CATEGORY_B=[]
CATEGORY_C=[]
CATEGORY_D=[]
CATEGORY_E=[]
CATEGORY_F=[]
def toString(numberString,courseName,list):

    numberString=numberString.split("; ")
    courses= [courseName + " "+x for x in numberString]
    courses=" ,".join(courses)
    print(courses)
    list.append(courses)
toString("1; 2; 3; 5; 6; 11; 15; 17; 18; 21; 22; 52; 71; 72","HIST",CATEGORY_A)
toString("2; 5","DANCE",CATEGORY_A)
toString("11; 45; 55","ENGL",CATEGORY_A)
toString("1; 2; 5; 6; 7; 8; 9; 11","FILM",CATEGORY_A)
toString("30; 31; 32; 33; 35; 36; 37; 39","MUSIC",CATEGORY_A)
toString("11","PHILOS",CATEGORY_A)
toString("52","PHOTO",CATEGORY_A)
toString("2; 5","TH ART",CATEGORY_A)
toString("3; 4; 5; 6; 7; 8; 9; 10; 14; 15; 17; 26; 34; 38; 39; 40; 41; 49; 50; 51; 52; 53; 54; 56; 57; 58; 59","ENGL",CATEGORY_B)
toString("9","CHNESE",CATEGORY_B)
toString("20","ENVRN",CATEGORY_B)
toString("1; 2; 3; 4; 5; 6; 11; 12; 19; 21; 22; 24; 25; 26; 29; 33; 34; 38; 39; 41; 45; 46; 48; 53; 55","HIST",CATEGORY_B)
toString("26","HUM",CATEGORY_B)
toString("9","JAPAN",CATEGORY_B)
toString("1; 2; 3; 4; 5; 6; 10; 20; 22; 23; 24; 48; 51; 52","PHILOS",CATEGORY_B)
toString("51; 52","POL SC",CATEGORY_B)
toString("51; 52","REL ST",CATEGORY_B)
toString("2; 3; 7; 14; 20; 21; 22","ANTHRO",CATEGORY_C)
toString("5; 6; 15","ECON",CATEGORY_C)
toString("7; 14; 22; 32","ENVRN",CATEGORY_C)
toString("2; 7; 8; 11; 14","GEOG",CATEGORY_C)
toString("5; 10; 11","GLOBAL",CATEGORY_C)
toString("10; 13; 14; 15; 16; 20; 28; 32; 42; 43; 52; 62","HIST",CATEGORY_C)
toString("3; 5; 7; 8; 14; 21; 22; 23; 31; 47","POL SC",CATEGORY_C)
toString("1; 1S; 2; 2S; 12; 30; 31; 32; 33; 34","SOCIOL",CATEGORY_C)
toString("8","URBAN",CATEGORY_C)
toString("10; 20; 30","WOM ST",CATEGORY_C)
toString("5","ANTHRO",CATEGORY_D)
toString("3; 4; 15; 21","BIOL",CATEGORY_D)
toString("3","PHYS",CATEGORY_D)
# toString("","",CATEGORY_D)
# toString("","",CATEGORY_D)
# toString("","",CATEGORY_D)
# toString("","",CATEGORY_D)
# toString("","",CATEGORY_D)
toString("3; 4","ASTRON",CATEGORY_E)
toString("9; 11; 19","CHEM",CATEGORY_E)
toString("5","GEOG",CATEGORY_E)
toString("4; 5","GEOL",CATEGORY_E)
toString("6; 8; 14; 21","PHYSCS",CATEGORY_E)
toString("1; 2","ECON",CATEGORY_F)
toString("2; 7; 21; 26; 28; 29; 54","MATH",CATEGORY_F)
toString("9","PHILOS",CATEGORY_F)
# toString("","",CATEGORY_E)
# toString("","",CATEGORY_F)
# toString("","",CATEGORY_F)
# toString("","",CATEGORY_F)
# toString("","",CATEGORY_F)
# toString("","",CATEGORY_F)
# toString("","",CATEGORY_F)
# toString("","",CATEGORY_F)
# toString("","",CATEGORY_F)
# toString("","",CATEGORY_F)

print(CATEGORY_A)
def getCategoryString(category):
    print(" ,".join(category)+" ")
    return " ,".join(category)+" "
getCategoryString(CATEGORY_A)
getCategoryString(CATEGORY_B)
getCategoryString(CATEGORY_C)
getCategoryString(CATEGORY_D)
getCategoryString(CATEGORY_E)
getCategoryString(CATEGORY_F)
print(getCategoryString(CATEGORY_A)+","+
getCategoryString(CATEGORY_B)+","+
getCategoryString(CATEGORY_C)+","+
getCategoryString(CATEGORY_D)+","+
getCategoryString(CATEGORY_E)+","+
getCategoryString(CATEGORY_F))