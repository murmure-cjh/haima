#!/haplox/users/donglf/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# reference: https://www.runoob.com/python/python-email.html
# reference: https://www.liaoxuefeng.com/wiki/1016959663602400/1017790702398272
# reference: https://docs.python.org/zh-cn/3/library/smtplib.html
# edit on 2024-07-25 by WX
# edit on 2025-03-12 by WX, add Wu LiuYu to sample_receviers
# edit on 2025-07-08 by WX, remove members from HaploX Shenzhen lab, add members from Jiangxi lab


import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.utils import parseaddr, formataddr
from email.header import Header
import os

mail_host="smtp.exmail.qq.com"
mail_user="chenjh@haplox.com"
#mail_pass="8XgtB92DKMU9R5Jm"
# mail_pass updated on 2024-07-25
# mail_pass="zr96cTRpTWEHEiK4"
# mail_pass updated on 2025-10-31
mail_pass="LiAdvMqR3vxujmkX"
mail_port=465
sender_name = 'chenjh'
sender_email = 'chenjh@haplox.com'

maintainer_receviers = ['chenjh@haplox.com','chenbf@haplox.com', 'linyushan@haplox.com', 'huyw@haplox.com', 'wangyz@haplox.com', 'wuliuyu@haplox.com']
maintainer_receviers_SZ = ['chenjh@haplox.com', 'linyushan@haplox.com', 'wangyz@haplox.com']
maintainer_receviers_JX = ['chenjh@haplox.com', 'wangyz@haplox.com', 'chenbf@haplox.com', 'huyw@haplox.com', 'wuliuyu@haplox.com']

start_email_subject = '【海码生产】批次 {} 启动'
start_email_content = '批次 {} 启动\n 详细信息请查看附件'

tcr_receviers = ['liangshu@haplox.com','gancw@haplox.com']
tcr_email_subject = '【海码TCR生产】批次 {} TCR样本信息'
tcr_email_content = '批次 {} TCR样本信息\n 详细信息请查看附件'

invalid_sample_receviers = ['wangyz@haplox.com', 'huyw@haplox.com', 'jxhgclab@haplox.com', 'wuliuyu@haplox.com', 'fanxh@haplox.com']
invalid_sample_receviers_SZ = ['wangyz@haplox.com', 'huyw@haplox.com', 'jxhgclab@haplox.com', 'wuliuyu@haplox.com', 'fanxh@haplox.com']
invalid_sample_receviers_JX = ['wangyz@haplox.com', 'huyw@haplox.com', 'jxhgclab@haplox.com', 'wuliuyu@haplox.com', 'fanxh@haplox.com']

invalid_email_subject = '【海码生产】批次 {} 不合格样本信息'
invalid_email_content = '批次 {} 不合格样本信息\n 详细信息请查看附件'

reporter_receviers = ['linyushan@haplox.com', 'chenbf@haplox.com']
reporter_receviers_SZ = ['linyushan@haplox.com']
reporter_receviers_JX = ['chenbf@haplox.com']
end_email_subject = '【海码生产】批次 {} 结束'
end_email_content = '批次 {} 结束\n 详细信息请查看附件'

test_receviers = ['chenjh@haplox.com']


def send_email(receviers, subject, content, attachment_file=None):
    """
    send email from sender to receviers with content and attachment files
    """
    try:
        message = MIMEMultipart()
        message['From'] = Header(f"{sender_name}<{sender_email}>", 'utf-8')
        message['To'] =  Header('; '.join(receviers), 'utf-8')
        message['Subject'] = Header(subject, 'utf-8')
        message.attach(MIMEText(f'{content}', 'plain', 'utf-8'))
        
        if attachment_file and os.path.exists(attachment_file):
            att_file = MIMEText(open(attachment_file, 'rb').read(), 'base64', 'utf-8')
            att_file_name = os.path.basename(attachment_file)
            att_file["Content-Type"] = 'application/octet-stream'
            att_file["Content-Disposition"] = f'attachment; filename={att_file_name}'
            message.attach(att_file)
            
        server = smtplib.SMTP_SSL(mail_host, mail_port)
        server.set_debuglevel(1)  # 开启详细调试信息
        server.login(mail_user, mail_pass)
        server.sendmail(sender_email, receviers, message.as_string())
        server.quit()
        print(f"Send mail successfully to {'; '.join(receviers)}")
    except Exception as e:
        print(f"Could not send mail to {'; '.join(receviers)}")
        print(f"Error type: {type(e).__name__}")
        print(f"Error details: {str(e)}")
        import traceback
        traceback.print_exc()  # 打印完整的错误堆栈

def send_email2(receviers, subject, content, attachment_files=None):
    """
    send email from sender to receviers with content and attachment files, multiple attachment_files
    """
    try:
        message = MIMEMultipart()
        message['From'] = Header(f"{sender_name}<{sender_email}>", 'utf-8')
        message['To'] =  Header('; '.join(receviers), 'utf-8')
        message['Subject'] = Header(subject, 'utf-8')
        message.attach(MIMEText(f'{content}', 'plain', 'utf-8'))
        for attachment_file in attachment_files:
            if not attachment_file:
                continue
            att_file = MIMEText(open(attachment_file, 'rb').read(), 'base64', 'utf-8')
            att_file_name = os.path.basename(attachment_file)
            att_file["Content-Type"] = 'application/octet-stream'
            att_file["Content-Disposition"] = f'attachment; filename={att_file_name}'
            message.attach(att_file)
        server = smtplib.SMTP_SSL(mail_host, mail_port)
        # server.set_debuglevel(1)
        server.login(mail_user, mail_pass)
        server.sendmail(sender_email, receviers, message.as_string())
        server.quit()
        print (f"Send mail successefully to {'; '.join(receviers)}")
    except:
        print (f"Could not send mail to {'; '.join(receviers)}")
        
        

if __name__ == "__main__":
    # 测试发送邮件
    send_email(
        receviers=['linyushan@haplox.com'], 
        subject='测试邮件',
        content='这是一封测试邮件',
        attachment_file=None 
    )