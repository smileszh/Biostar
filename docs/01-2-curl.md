


# curl
curl 是一个命令行工具和库，用于传输数据，支持多种协议，包括 HTTP、HTTPS、FTP、FTPS、SFTP、IMAP、SMTP、POP3、LDAP、RTMP 和 RTSP。curl 还支持 SSL 证书、HTTP POST、HTTP PUT、FTP 上传、HTTP 基本身份验证、代理、cookie、用户代理、压缩、断点续传、文件传输限速、重定向等功能。

使用示范:

``` bash
URL=https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

make -f src/run/curl.mk URL=${URL} run
```


[aria.mk](src/run/curl.mk) 


