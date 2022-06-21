# SoftGNSS
Current working ver. of SoftGNSS v3.0 for GN3sV2, GN3sV3, NT1065EVK, and NUT4NT samplers. All known updates included. Navigation module updated.

## Note 
This code has been adapted by Xin Zhang for purposes of course "AV423 Satellite Navigation" taught at School of Aeronautics & Astronautics, Shanghai Jiao Tong University, from the SoftGNSS v3.0 code base developed for the text: "A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach" by Borre, Akos, et.al.
This repo and my entire research will not be possible without Prof Borre's (I am sorry that he has passed) EASY suite, the seminal work on GPS software receiver by Prof Akos, and their generous willingness to share SoftGNSS v3.0. 
Please check Google Scholar for a complete list of Prof. Akos's recent publications and patents. 



## Usage
"Straight out of the box": runs with GN3sV2

or using data file downloaded here:

gnss0.bin:
[Baidu Netdisk](https://pan.baidu.com/s/11dcgJer-Cz4HcFvoWQ_8_w) PIN: crvi

or [Googledrive](https://drive.google.com/file/d/1QJzBehlsqkjN_NX6GtBDtODARB5SS1Yb/view?usp=sharing)

gnsa14.bin:
[Baidu Netdisk](https://pan.baidu.com/s/1Kla89mZ4XCIzBcNKpKVYng). PIN: 61p0

or [GoogleDrive](https://drive.google.com/file/d/1SEveFaYKCw9DlSlD-XE5Fc3QWeILJ45T/view?usp=sharing)

The two data files were collected using [SiGe GN3S Sampler v2](https://www.sparkfun.com/products/retired/8238) (OOP though...) and therefore,

|File Name   | Data Length in Seconds | I\/Q Format (as in initSettings.m)        |
| ---------- | ---------------------- | ----------------------------------------- |
| gnss0.bin  | 60.5                   | 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,I3/Q3 |
| gnsa14.bin | 40                     | 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,I3/Q3 |


## Update
This receiver will receive a major update in one month (as of 2022-6-22 UTC+8...oops sorry...).
