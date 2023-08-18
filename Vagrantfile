# -*- mode: ruby -*-
# vi: set ft=ruby :

required_plugins = %w( vagrant-vbguest )
required_plugins.each do |plugin|
  system "vagrant plugin install #{plugin}" unless Vagrant.has_plugin? plugin
end

Vagrant.configure("2") do |config|
  # You can search for boxes at https://app.vagrantup.com/fedora
  config.vm.box = "fedora/34-cloud-base"
  config.ssh.forward_agent = true
  config.ssh.forward_x11 = true

  config.vm.synced_folder ".", "/vagrant", disabled: true
  config.vm.synced_folder ".", "/home/vagrant/flow-simulator", rsync__args: ["--rsync-path='sudo rsync'", "--archive"]
  config.vm.network "forwarded_port", guest: 9080, host: 9080

  config.vm.provider "virtualbox" do |vb|
    vb.name = "flow-simulator"
    vb.gui = true
    vb.memory = "4096"
    vb.cpus = "4"
  end

  config.vm.provider "libvirt" do |lv|
    lv.memory = "4096"
    lv.cpus = "4"
  end

  config.vm.provision "shell", inline: <<-SHELL
    # Setting up locale configurations
    echo "LANG=en_US.UTF-8"                         >> /etc/environment
    echo "LANGUAGE=en_US.UTF-8"                     >> /etc/environment
    echo "LC_ALL=en_US.UTF-8"                       >> /etc/environment
    echo "LC_CTYPE=en_US.UTF-8"                     >> /etc/environment
    echo "PYTHONDONTWRITEBYTECODE=1"                >> /etc/environment
    echo "ROOT_DIR=/home/vagrant/flow-simulator"    >> /etc/environment

    # Setting up bash environment
    echo "alias cmake-generate='cmake -GNinja -Bbuild -DCMAKE_BUILD_TYPE=Release -H.'"         >> /etc/profile
    echo "alias cmake-generate-profile='cmake -GNinja -Bbuild -DCMAKE_BUILD_TYPE=Profile -H.'" >> /etc/profile
    echo "alias cmake-build='cmake --build build/ --target flow-simulator.x'"                       >> /etc/profile
    echo "alias cmake-build-test='cmake --build build/ --target flow-simulator-test.x'"             >> /etc/profile
    echo "alias cmake-build-profile='cmake --build build/ --target flow-simulator-profile.x'"       >> /etc/profile
    echo "alias cmake-build-all='cmake --build build/ --target all'"                           >> /etc/profile
    echo "alias cmake-clean='cmake --build build/ --target clean'"                             >> /etc/profile
    echo "alias cmake-clean-all='cmake --build build/ --target clean-all'"                     >> /etc/profile
    echo "cd ${ROOT_DIR}"                                                                      >> /etc/profile

    # Install flow-simulator dependencies
    echo "Installing dependencies..."
    dnf --best --allowerasing -y install \
      armadillo-devel \
      cmake \
      findutils \
      gcc-c++ \
      glog-devel \
      gperftools \
      gtest-devel \
      lcov \
      libffi-devel \
      ninja-build \
      pipenv \
      pkgconfig \
      rapidjson \
      redhat-rpm-config \
      redis \
      sundials-devel \
      tclap \
      && dnf clean all

    echo "\nEnter the command below to install Python dependencies..."
    echo "$ sudo pipenv install --system --skip-lock --clear"
  SHELL
end
