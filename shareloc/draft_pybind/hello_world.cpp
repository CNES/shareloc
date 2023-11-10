#include <pybind11/pybind11.h>
#include "hello_world.h"

HW::HW() {}

HW::~HW() {}

std::string HW::hellow_world() const {
    return "Hello world !";
}

/*int main(){
    HW HelloWorld;
    std::cout<<HelloWorld.hellow_world()<<std::endl;
}*/