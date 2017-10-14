# MSstatsTMT

Before build/install, remember to updata .rd file by useing "devtools::document()", also in Test.R file (Make sure you are in the package directory.

#Build package

Local build and install
1. Go to directory where the package fold located
2. Type: R CMD build MSstatsTMT (normally this will generate a file with name like "MSstatsTMT_0.0.0.9000.tar.gz"
3. Type: R CMD INSTALL MSstatsTMT_0.0.0.9000.tar.gz (or what the package called)


Github install

install_github("Vitek-Lab/MSstatsTMT/MSstatsTMT")
(currently not working due to private repo)

