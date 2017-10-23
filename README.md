# MSstatsTMT

Before build/install, remember to updata .rd file by useing "devtools::document()", also in Test.R file (Make sure you are in the package directory.

#Build package

Local build and install
1. Go to directory where the package folder located
2. Type: R CMD build MSstatsTMT (normally this will generate a file with name like "MSstatsTMT_0.0.0.9000.tar.gz"
3. Type: R CMD INSTALL MSstatsTMT_0.0.0.9000.tar.gz (or whatever the package is called)


Github install

devtools::install_github(repo = "Vitek-Lab/MSstatsTMT",auth_token = "0c4ece2b7e9aab30c8b1941e319a2c97fb2eb94e")


Note: Please let Sicheng know if the token is not working, and he will fix it. Also, don't share the token.

