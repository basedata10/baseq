#!/usr/bin/env python
# -*- coding: utf-8 -*-

import smtplib, os
from email.mime.multipart import MIMEMultipart
from email.mime.image import MIMEImage
from email.mime.text import MIMEText
from email.header import Header
from baseq.mgt import config

class Client:
    def __init__(self):
        # 邮件信息
        self.server = 'smtp.163.com'
        self.sender = 'baseq_admin@163.com'
        #self.receiver = ['1136177526@qq.com', '1606390004@pku.edu.cn', 'friedpine@gmail.com']
        self.receiver = ['1606390004@pku.edu.cn']
        self.password = 'shouquanqwer1234'

    def send_mail(self, subject="BaSeq Results", msg="", attches=[]):
        # 设置邮件内容
        msgRoot = MIMEMultipart('mixed')
        msgRoot['Subject'] = subject
        msgRoot['From'] = Header("BeiSeq Server", 'utf-8')
        msgRoot['To'] = ','.join(self.receiver)

        # 构造正文
        mail_msg = """  
            <h2>这是标题2</h2>
            <h3>这是标题3。。。</h3>              
            <p>Python 邮件发送测试...</p>  
            <p><a href="http://www.runoob.com">这是一个链接</a></p>
            """
        message = MIMEText(mail_msg, 'html', 'utf-8')
        msgRoot.attach(message)

        # 构造附件
        for file in attches:
            if file[-3:] == "png":
                 att = MIMEImage(open(file, 'rb').read())
                 att.add_header('Content-ID', '<{}>'.format(file))
                 att.add_header("Content-Type", "image/png")
                 att.add_header('Content-Disposition', 'attachment; filename="'+file[0:-4]+'"')
            else:
                att = MIMEText(open(file, 'rb').read(), 'base64', 'utf-8')
                att["Content-Type"] = 'application/octet-stream'
                att["Content-Disposition"] = 'attachment; filename="'+file+'"'
            msgRoot.attach(att)

        # 发送邮件
        try:
            print("[info] Start Connect SMTP")
            smtp = smtplib.SMTP(self.server)
            smtp.set_debuglevel(1)
            print("[info] Start Logging In...")
            print(self.sender, self.password)
            smtp.login(self.sender, self.password)
            print("[info] Success Logging In...")
            smtp.sendmail(self.sender, self.receiver, msgRoot.as_string())
            smtp.quit()
            print("send mail succeed")
        except:
            print("send mail error")

# 测试
if __name__ == '__main__':
    msg = ""
    client = Client()
    #client.send_mail(subject, message, attches=[])