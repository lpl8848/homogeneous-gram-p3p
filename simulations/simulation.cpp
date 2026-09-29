#include <gram_p3p/p3p.h>
#include "data.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace fs=std::filesystem;
using Pose=gram_p3p::Pose;
using Scene=simulation::Scene;

double error(const Pose& p,const Scene& s) {
    return (p.R()-s.R).cwiseAbs().sum()+(p.t-s.t).cwiseAbs().sum();
}
bool valid(const Pose& p,const Scene& s,bool sphere) {
    if (!p.R().allFinite() || !p.t.allFinite()) return false;
    if (std::abs(p.R().determinant()-1)>1e-6 ||
        (p.R().transpose()*p.R()-Eigen::Matrix3d::Identity()).cwiseAbs().sum()>1e-6 ||
        std::abs(Eigen::Quaterniond(p.R()).norm()-1)>1e-5) return false;
    double residual=0;
    for (int i=0;i<3;++i) {
        const Eigen::Vector3d point=p.R()*s.world[i]+p.t;
        if (point.dot(s.bearings[i])<=0) return false;
        if (sphere) {
            const double scale=point.dot(s.bearings[i]);
            residual+=(point-scale*s.bearings[i]).norm()/std::max(point.norm(),1e-300);
        } else {
            if (point.z()<=0) return false;
            residual+=(point.hnormalized()-s.bearings[i].hnormalized()).cwiseAbs().sum();
        }
    }
    return std::isfinite(residual) && residual <= (sphere?1e-7:1e-4);
}

void run(const fs::path& directory,std::uint64_t samples,std::uint64_t seed,
         double epsilon,bool sphere,std::uint64_t save_samples) {
    if (fs::exists(directory)) throw std::runtime_error("Output directory exists; choose a fresh path.");
    fs::create_directories(directory);
    std::ofstream errors(directory/"errors.csv"), failures(directory/"failure_inputs.csv"),
                  counts(directory/"summary.csv"), config(directory/"config.json");
    if (!errors || !failures || !counts || !config) throw std::runtime_error("Cannot write results.");
    errors<<std::setprecision(17)<<"id,best_pose_L1,returned,correct,unique_correct\n";
    failures<<std::setprecision(17)<<"id,best_pose_L1";
    for (const char* prefix : {"f","X","R"})
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) failures<<','<<prefix<<i<<j;
    failures<<",t0,t1,t2\n";
    simulation::Generator generator(seed);
    std::vector<Pose> poses;poses.reserve(4);
    std::uint64_t returned=0,correct=0,unique=0,duplicates=0,incorrect=0,
        gt=0,empty=0,no_correct=0,failure_records=0;
    for (std::uint64_t i=0;i<samples;++i) {
        const Scene scene=generator.draw(epsilon,sphere);
        gram_p3p::solve(scene.bearings,scene.world,&poses);
        double best=std::numeric_limits<double>::infinity();
        std::vector<std::size_t> distinct;int good=0;
        for(std::size_t j=0;j<poses.size();++j) {
            const Pose& pose=poses[j];const double e=error(pose,scene);
            if (std::isfinite(e)) best=std::min(best,e);
            if (!valid(pose,scene,sphere)) {++incorrect;continue;}
            ++good;bool duplicate=false;
            for (auto k:distinct) if ((pose.R()-poses[k].R()).cwiseAbs().sum()+
                    (pose.t-poses[k].t).cwiseAbs().sum()<1e-5) duplicate=true;
            if (!duplicate) distinct.push_back(j);
        }
        returned+=poses.size();correct+=good;unique+=distinct.size();
        duplicates+=good-distinct.size();empty+=poses.empty();no_correct+=distinct.empty();
        gt+=best<=1e-6;
        if(i<save_samples) errors<<i<<','<<best<<','<<poses.size()<<','<<good<<','<<distinct.size()<<'\n';
        // A bounded file size still retains the first 10,000 failure inputs.
        if(best>1e-6 && failure_records<10000) {
            failures<<i<<','<<best;
            for(const auto& v:scene.bearings) for(int j=0;j<3;++j) failures<<','<<v[j];
            for(const auto& v:scene.world) for(int j=0;j<3;++j) failures<<','<<v[j];
            for(int a=0;a<3;++a) for(int b=0;b<3;++b) failures<<','<<scene.R(a,b);
            for(int j=0;j<3;++j) failures<<','<<scene.t[j];failures<<'\n';++failure_records;
        }
    }
    counts<<"samples,seed,epsilon,returned,correct,unique_correct,duplicates,incorrect,GT_recovered,GT_failures,empty,no_correct,failure_records\n"
          <<std::setprecision(17)<<samples<<','<<seed<<','<<epsilon<<','<<returned<<','<<correct<<','<<unique
          <<','<<duplicates<<','<<incorrect<<','<<gt<<','<<samples-gt<<','<<empty<<','<<no_correct<<','<<failure_records<<'\n';
    config<<std::setprecision(17)<<"{\n  \"samples\": "<<samples<<",\n  \"seed\": "<<seed
          <<",\n  \"epsilon\": "<<epsilon<<",\n  \"full_sphere\": "<<(sphere?"true":"false")
          <<",\n  \"saved_error_samples\": "<<std::min(samples,save_samples)
          <<",\n  \"max_failure_records\": 10000,\n  \"GT_threshold\": 1e-6\n}\n";
    std::cout<<directory.string()<<": GT recovered "<<gt<<'/'<<samples
             <<", incorrect outputs "<<incorrect<<'\n';
}

int main(int argc,char** argv) {
    try {
        std::string mode="accuracy";fs::path output="results/accuracy";
        std::uint64_t samples=100000,seed=1,save_samples=100000;
        for(int i=1;i<argc;++i) {
            const std::string option=argv[i];
            if(option=="--help") {std::cout<<"--mode accuracy|collinear|sphere --samples N --seed S --output DIR --save-samples N\n";return 0;}
            if(i+1==argc) throw std::runtime_error("Missing option value.");
            const std::string value=argv[++i];
            if(option=="--mode") mode=value;
            else if(option=="--output") output=value;
            else if(option=="--samples") samples=std::stoull(value);
            else if(option=="--seed") seed=std::stoull(value);
            else if(option=="--save-samples") save_samples=std::stoull(value);
            else throw std::runtime_error("Unknown option: "+option);
        }
        if(samples==0 || (mode!="accuracy" && mode!="collinear" && mode!="sphere"))
            throw std::runtime_error("Use a valid mode and positive sample count.");
        if(mode=="collinear") {
            if(fs::exists(output)) throw std::runtime_error("Output directory exists; choose a fresh path.");
            for(int power=2;power<=10;++power)
                run(output/("e"+std::to_string(power)),samples,seed+power,
                    std::pow(10.,-power),false,save_samples);
        } else run(output,samples,seed,0,mode=="sphere",save_samples);
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
    return 0;
}
