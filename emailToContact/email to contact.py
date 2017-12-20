import csv
emailList=[]
with open ("/media/toor/storage/web/script/emailToContact/email.csv", newline="") as csvfile:
    reader = csv.reader(csvfile,delimiter=",")
    for x in reader:
        emailList.append(",,,,,,,,,,,,,,,,,,,,,,,,,,,* ,"+x[0])

with open ("/media/toor/storage/web/script/emailToContact/contact.csv","w",newline="") as contactFile:
    # spamwriter = csv.writer(contactFile, delimiter='\t')
    contactFile.write("Name,Given Name,Additional Name,Family Name,Yomi Name,Given Name Yomi,Additional Name Yomi,Family Name Yomi,Name Prefix,Name Suffix,Initials,Nickname,Short Name,Maiden Name,Birthday,Gender,Location,Billing Information,Directory Server,Mileage,Occupation,Hobby,Sensitivity,Priority,Subject,Notes,Group Membership,E-mail 1 - Type,E-mail 1 - Value\n")
    for x in emailList:
        contactFile.write(x+"\n")


