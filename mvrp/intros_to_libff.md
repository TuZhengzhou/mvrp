libff 库介绍

  - common 文件夹：

    - (1) 采用计时

      libff/common/profiling.hpp

      libff::start_profiling();

    - (2) 
  
  - algebra/curves 文件夹

    - (1) 公共参数

      libff/algebra/curves/public_params.hpp

      libff::G1<ppT>
      
      libff::G2<ppT>

      libff::GT<ppT>

      libff::Fr<ppT>

    - (2) mnt 曲线

      libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp

      使用首先需要初始化参数：libff::mnt6_pp::init_public_params()

    - (3) bn_128 曲线

      libff/algebra/curves/bn128/bn128_pp.hpp

      使用首先需要初始化参数：
      libff::bn128_pp::init_public_params()

      