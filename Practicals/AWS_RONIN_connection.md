# Connecting to your Coud Compute Environment

You will have each been given access to cloud compute resources which we access using the IP address below.
All cloud computing accounts are private but are identical and we will use these for the entire practical series. 
Accessing your cloud compute resource is like having your very own server. 

**Follow all instructions here very, very carefully**

1. Open your favourite internet browser. We recommend Firefox, but Edge/Chrome are also acceptable. Safari has not been tested.
2. Enter the following address (`https://rstudio-ubuntu.uoa.cloud:4200`) in the address bar of your browser.
3. You should see this login screen  

    ![AWS RONIN shell in a box](./Bash_Practicals/images/shell_in_a_box.png)

4. Login with your UoA id number and password.

    Now you should see the following: 

    ![Rstudio_pwd_from_shell](./Bash_Practicals/images/RStudio_server_spawn.png)

    Once the Rstudio server has been launched (this may take a couple of minutes) you should see the following message in the browser: 

    ```
    Welcome to rstudio-ubuntu.uoa.cloud

    You can now access Rstudio by opening a new browser and type

    http://rstudio-ubuntu.uoa.cloud:8001

    Here are the credentials for this session : !

    username: a1234567 (your id here)

    password: (a 20 character password for your Rstudio login)
    ```
    __Make sure you copy the password before proceeding__

5. Open a new browser tab as instructed above and type in the URL given above (`http://rstudio-ubuntu.uoa.cloud:8001`)*For this step it is essential that you use `http` and NOT `https`*

    you should see the following login screen. 

    ![Rstudio_login_screen](./Bash_Practicals/images/Rstudio_AWS_login.png)

6. Use your UoA id (a1234567) and paste the password that you copied in the password box. __You have two minutes to complete this__ if you don't complete in time, you have to start over. 

    You now have access to your cloud compute resource (VM). 

    - Every time you log in to your VM you will be given a new one time password that you have to paste into the Rstudio login panel.
    - If you leave Rstudio idle for too long the VM will shut down. This is because we pay by the minute for access to this resource. As long as you have jobs running or do not leave the terminal idle for more than an hour your VM will not automatically shut down. 

    Importantly, because of the onetime password you will be the only person to be able to access your VM as the system has been set up to allow you to use one VM at a time.

You can access your VM:

- When connected directly to the University WiFi Network
- When connected from another network using the [University of Adelaide VPN](https://www.adelaide.edu.au/technology/your-services/network-services/remote-access-via-virtual-private-network-vpn)

[Back](./Bash_Practicals/1_IntroBash.md)