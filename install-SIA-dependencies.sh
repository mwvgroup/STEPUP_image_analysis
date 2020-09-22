#!/bin/bash


# Make sure python3 is installed
if ! command -v python3 &> /dev/null
then
    echo -e "\e[1;31minstallation failed:  python3 could not be found."
    echo -e "Make sure python3 is installed and in your PATH.\e[0m"
    exit
fi

# Make sure pip3 is installed
if ! command -v pip3 &> /dev/null
then
    echo -e "\e[1;31minstallation failed:  pip3 could not be found."
    echo -e "Make sure pip3 is installed.\e[0m"
    exit
fi


# Pip install the necessary python python3 packages
pip3 install matplotlib numpy astropy photutils scipy


# Install wcstools 3.9.6 if not already installed
if [ ! -d "/usr/local/wcstools-3.9.6" ]
then
	# Download wcstools 3.9.6 from the internet
	sudo curl  http://tdc-www.harvard.edu/software/wcstools/wcstools-3.9.6.tar.gz -o /usr/local/wcstools-3.9.6.tar.gz

	cd /usr/local/

	# Expand the download
	sudo tar -xf wcstools-3.9.6.tar.gz -C /usr/local/

	# Remove the archived version of wcstools
	sudo rm wcstools-3.9.6.tar.gz

	cd wcstools-3.9.6/

	# Compile the wcstools stuff
	sudo make all
fi


# Add wcstools to PATH if not already there
if ! grep -c "export PATH=\$PATH:/usr/local/wcstools-3.9.6/bin" ~/.bashrc
then
        echo "export PATH=\$PATH:/usr/local/wcstools-3.9.6/bin" >> ~/.bashrc
fi

if ! grep -c "export PATH=\$PATH:/usr/local/wcstools-3.9.6/bin" ~/.zshrc
then
        echo "export PATH=\$PATH:/usr/local/wcstools-3.9.6/bin" >> ~/.zshrc
fi


echo -e "\e[1;31m"
echo "RESTART YOUR TERMINAL/SHELL BEFORE RUNNING SIA!"
echo "If you do not restart your terminal before running, it will not work"
echo -e "\e[0m"
