# -*- coding: utf-8 -*-
"""
Created on Tue May  9 23:48:27 2017

@author: Hughie
"""

#***************获得男女比例*******************************
import itchat
f = open('signature.txt','w')
itchat.login()#登录
friends = itchat.get_friends(update=True)[0:]# 获取好友列表
male = female = other = 0# 初始化计数器，有男有女，当然，有些人是不填的
for i in friends[1:]:# 遍历这个列表，列表里第一位是自己，所以从"自己"之后开始计算# 1表示男性，2女性
    sex = i["Sex"]
    if sex == 1:
        male += 1
    elif sex == 2:
        female += 1
    else:
        other += 1# 总数算上，好计算比例啊～
total = len(friends[1:])
print('You have '+str(total)+' friends')
print (u"Male friends：%.2f%%" % (float(male) / total * 100))
print (u"Female friends：%.2f%%" % (float(female) / total * 100))
print (u"Others：%.2f%%" % (float(other) / total * 100))

#***************图表化男女比例*******************************
'''
from echarts import Echart,Legend, Piechart = Echart(u'%s的微信好友性别比例' % (friends[0]['NickName']), 'from WeChat')
chart.use(Pie('WeChat',[{'value': male, 'name': u'男性 %.2f%%' % (float(male) / total * 100)},{'value': female, 'name': u'女性 %.2f%%' % (float(female) / total * 100)},{'value': other, 'name': u'其他 %.2f%%' % (float(other) / total * 100)}],radius=["50%", "70%"]))
chart.use(Legend(["male", "female", "other"]))
del chart.json["xAxis"]
del chart.json["yAxis"]
chart.plot()
'''
#***************好友签名列表获取*******************************
'''
friends = itchat.get_friends(update=True)[0:]
for i in friends:# 获取个性签名
    signature = i["Signature"].strip().replace("span", "").replace("class", "").replace("emoji", "")# 正则匹配过滤掉emoji表情，例如emoji1f3c3等
    rep = re.compile("1f\d.+")
    signature = rep.sub("", signature)
    print(str(signature))
'''
#***************自动回复消息*******************************
# 封装好的装饰器，当接收到的消息是Text，即文字消息
@itchat.msg_register('Text')

def text_reply(msg):
    # 当消息不是由自己发出的时候
    if not msg['FromUserName'] == myUserName:
        # 发送一条提示给文件助手
        itchat.send_msg(u"[%s]收到好友@%s 的信息：%s\n" %
                        (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(msg['CreateTime'])),
                         msg['User']['NickName'],
                         msg['Text']), 'filehelper')
        # 回复给好友
        return u'[自动回复]您好，我现在有事不在，一会再和您联系。\n已经收到您的的信息：%s\n' % (msg['Text'])
if __name__ == '__main__':
    itchat.auto_login()
    myUserName = itchat.get_friends(update=True)[0]["UserName"] # 获取自己的UserName
    itchat.run()        