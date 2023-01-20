/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */

#ifndef CODECFACTORY_H_
#define CODECFACTORY_H_

#include "common.h"
#include "codecs.h"
#include "util.h"
#include "FOR.h"
#include "FOR_my.h"
#include "combinedcodec.h"
#include "piecewise.h"
#include "piecewise_fiting.h"
#include "piecewise_fix.h"
#include "piecewise_cost.h"
#include "piecewise_cost_float.h"
#include "piecewise_cost_lookahead.h"
#include "piecewise_cost_dp.h"
#include "piecewise_fix_op.h"
#include "piecewise_fix_op_lp.h"
#include "piecewise_fix_op_lp_cost.h"
#include "piecewise_fix_op_minimize_maxerror.h"
#include "piecewise_fix_op_minimize_maxerror_round.h"
#include "piecewise_fix_op_round.h"
#include "piecewise_fix_op_float.h"
#include "piecewise_fix_merge.h"
#include "piecewise_fix_pack.h"
#include "piecewise_ransac.h"
#include "piecewise_double.h"
#include "piecewise_outlier_detect.h"
#include "ransac_outlier_detect.h"
#include "ransac_fix.h"
#include "piecewise_varilength.h"
#include "rle.h"
#include "nonlinear_fix.h"
#include "spline_fix.h"
#include "piecewise_fanout.h"
#include "piecewise_multi_fanout.h"
#include "maskvbyte.h"
#include "variablebyte.h"
#include "snappy.h"
#include "delta_my.h"
#include "piecewise_fix_delta.h"
#include "spline_fix.h"


namespace Codecset {

typedef std::map<std::string, IntegerCODEC*> CodecMap;

/**
 * This class is a convenience class to generate codecs quickly.
 * It cannot be used safely in a multithreaded context where
 * each thread should have a different IntegerCODEC.
 */
class CODECFactory {
public:
  static CodecMap scodecmap;

  // hacked for convenience
  static std::vector<IntegerCODEC*> allSchemes() {
    std::vector<IntegerCODEC*> ans;
    for (auto i = scodecmap.begin(); i != scodecmap.end(); ++i) {
      ans.push_back(i->second);
    }
    return ans;
  }

  static std::vector<std::string> allNames() {
    std::vector<std::string> ans;
    for (auto i = scodecmap.begin(); i != scodecmap.end(); ++i) {
      ans.push_back(i->first);
    }
    return ans;
  }

  static IntegerCODEC* &getFromName(std::string name) {
    if (scodecmap.find(name) == scodecmap.end()) {
      std::cerr << "name " << name << " does not refer to a CODEC."
                << std::endl;
      std::cerr << "possible choices:" << std::endl;
      for (auto i = scodecmap.begin(); i != scodecmap.end(); ++i) {
        std::cerr << static_cast<std::string>(i->first)
                  << std::endl; // useless cast, but just to be clear
      }
      std::cerr << "for now, I'm going to just return 'copy'" << std::endl;
      return scodecmap["copy"];
    }
    return scodecmap[name];
  }
};

// C++11 allows better than this, but neither Microsoft nor Intel support C++11
// fully.
static inline CodecMap initializefactory() {
  CodecMap map;

  map["FOR"]= new Codecset::FOR();
  map["FOR_my"]= new FOR_my();
  map["rle"]= new rle();
  map["delta_my"] = new delta_my();
  map["piecewise_fiting"]= new piecewise_fiting();
  map["piecewise_cost"]= new piecewiseCost();
  map["piecewise_cost_float"]= new piecewiseCost_float();
  map["piecewise_cost_dp"]= new piecewiseDp();
  map["piecewise_fix_op"]= new piecewise_fix_op();
  map["piecewise_fix_op_lp"]= new piecewise_fix_op_lp();
  map["piecewise_fix_op_lp_cost"]= new piecewise_fix_op_lp_cost();
  map["piecewise_fix_op_max"]= new piecewise_fix_op_max();
  map["piecewise_fix_op_max_round"]= new piecewise_fix_op_max_round();
  map["piecewise_fix_op_round"]= new piecewise_fix_op_round();
  map["piecewise_fix_op_float"]= new piecewise_fix_op_float();
  map["piecewise_fix_merge"]= new piecewise_fix_merge();
  map["MaskVByte"] = new Codecset::MaskVByte();
  map["spline_fix"] = new spline_fix();
  // map["snappy"] = new Codecset::JustSnappy();
  map["copy"] = new Codecset::JustCopy();

  return map;
}

CodecMap CODECFactory::scodecmap = initializefactory();

} // namespace FastPFor

#endif /* CODECFACTORY_H_ */
