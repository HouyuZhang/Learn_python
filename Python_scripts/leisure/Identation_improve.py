###################################################################
#This srcipt is for improving indentation for given C source codes
#Date:2016-12-16
#Version:1.0
###################################################################
import re
File=input('Please input the filename you want processing:\n')
f=open(File)      #Open the source code 
k=open('After_'+File,'w')#Write the after-processing source code to a new file

def rs():                   #This function is for reading and pre-processing
	line=f.readline()         #we can get a line from file
	p=re.compile('^[\t ]*')   #
	line=p.sub('',line)       #remove space and tab before the get line
	return line
	
def define_find(line):                                                        #This function is for finding #define without brace like this:
	if(re.match('#define.*\(.*\).*\(\s*',line) and not bool(re.search('{',line))): #   #define neq(X,Y) (
		k.write(line)                                                             #     .......
		line=rs()                                                                 #   )
		while(bool(re.match('\)',line))==False):
			k.write('\t'+line)
			line=rs()
		k.write(line)
		return 1

def decrease_tab(line):#This function is for finding '}' to reduce tabs
	global flag
	if(re.match('(}|(/\*.*\*/\s*}.*))',line) and not bool(re.match('(}\s*}.*)|(}.*{)',line))):   #only for finding this'}'   we can reduce one tab
		flag-=1
		k.write('\t'*flag+line)
		return 1
	if(re.match('}.*}.*',line)):                                                                 #only for finding this'}.*}' we can reduce two tabs
		flag-=2
		k.write('\t'*flag+line)
		return 1
	if(re.match('}.*{',line)):                                                                   #only for finding this'}.*{'  we can't reduce tab
		k.write('\t'*(flag-1)+line) 
		return 1

def increase_tab(line):#This function is for finding '{' or for(),while(),if() to add tabs
	global flag
	if(re.search('{',line) and not bool(re.search('({.*})|}',line))):                      #finding '{' but without '}' to add a tab
		k.write('\t'*flag+line)
		flag+=1
		return 1
	elif(re.search('(for|while|if)\s+\(.*\)\s*',line) and not bool(re.search('{|;',line))):#finding for|while|if(...) without '{'
		k.write('\t'*flag+line)
		flag+=1
		line=rs()
		if(re.match('{',line) and not bool(re.search(';',line))): #this small part if for a form like this:
			flag-=1                                                 #void func(...)
			k.write('\t'*flag+line)                                 #{
			flag+=1                                                 #	...
			return 1                                                #}
		elif(increase_tab(line)):
			pass
			return 1
		elif(re.search(';',line)):
			k.write('\t'*flag+line)
			flag-=1
		return 1
	
flag=0
line=rs()
while line:#For every given line we call these three functions respectively
	if(define_find(line)):
		pass
	elif(increase_tab(line)):
		pass
	elif(decrease_tab(line)):
		pass
	else:
		k.write('\t'*flag+line)
	line=rs()

f.close()
k.close()
