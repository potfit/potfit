/****************************************************************
 *
 * kim_read_model_props.cpp: list KIM model params for testing
 *
 ****************************************************************
 *
 * Copyright © 2018 the potfit development team
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the “Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall
 * be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY
 * KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
 * AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * https://www.potfit.net/
 *
 ***************************************************************/

#include <KIM_SimulatorHeaders.hpp>

#include <iostream>
#include <functional>
#include <map>
#include <memory>
#include <string>

namespace {
  constexpr auto g_usage = "Usage: kim_read_model_props <kim_model_name>";
  constexpr auto g_delimiter = '#';

  using kim_model_p = std::unique_ptr<KIM::Model, std::function<void(KIM::Model*)>>;

  std::ostream& operator<< (std::ostream& stream, const KIM::DataType& dt) {
    stream << "\"" << dt.ToString() << "\"";
    return stream;
  }

  template <class Head, class... Tail>
  void print1(const Head& head, Tail const&... tail){
    if constexpr (std::is_same<std::string, Head>::value)
      std::cout << "\"" << head << "\"";
    else
      std::cout << head;
    if constexpr (sizeof...(tail) > 0) {
      std::cout << g_delimiter;
      print1(tail...);
    }
  }

  template <class... Args, typename T>
  void print_item(T name, const Args&... args)
  {
    std::cout << name << "=";
    print1(args...);
    std::cout << "\n";
  }

  template <typename T, typename IT, typename std::enable_if<std::is_same<typename std::iterator_traits<IT>::iterator_category, std::random_access_iterator_tag>::value, IT>::type* = nullptr>
  void print_item(T name, IT begin, IT end)
  {
    std::cout << name << "=";
    for (auto val = begin; val != end; ++val) {
      if constexpr (std::is_same<std::string, typename std::iterator_traits<IT>::value_type>::value)
        std::cout << "\"" << *val << "\"";
      else
        std::cout << *val;
      if (val != (end - 1))
        std::cout << g_delimiter;
    }
    std::cout << "\n";
  }

  auto create_kim_model(const std::string& model_name)
  {
    KIM::Model* pmodel_temp = nullptr;

    int requestedUnitsAccepted {0};
    auto ret = KIM::Model::Create(KIM::NUMBERING::zeroBased, KIM::LENGTH_UNIT::A, KIM::ENERGY_UNIT::eV,
                                  KIM::CHARGE_UNIT::e, KIM::TEMPERATURE_UNIT::unused, KIM::TIME_UNIT::unused,
                                  model_name, &requestedUnitsAccepted, &pmodel_temp);

    if (ret)
      throw std::runtime_error("Error creating KIM_Model " + model_name);

    if (!requestedUnitsAccepted)
      throw std::runtime_error("KIM_Model " + model_name + " did not except standard params");

    print_item("NAME", model_name);

    return pmodel_temp;
  }

  void print_parameters(int num, kim_model_p& pmodel)
  {
    for (auto i = 0; i < num; ++i) {
      KIM::DataType dt;
      int extent = 0;
      const std::string* name;
      const std::string* desc;
      auto ret = pmodel->GetParameterMetadata(i, &dt, &extent, &name, &desc);
      if (ret)
        throw std::runtime_error("Error reading parameter " + std::to_string(i));
      print_item(std::string("PARAM_") + std::to_string(i), *name, dt, extent, *desc);
    }
  }

  void print_species(kim_model_p& pmodel)
  {
    int numberOfSpeciesNames = 0;

    KIM::SPECIES_NAME::GetNumberOfSpeciesNames(&numberOfSpeciesNames);

    std::vector<std::string> species;

    for (int i = 0; i < numberOfSpeciesNames; ++i) {
      KIM::SpeciesName sname;
      auto ret = KIM::SPECIES_NAME::GetSpeciesName(i, &sname);
      if (ret)
        throw std::runtime_error("Error reading species name for index " + std::to_string(i));
      int supported = 0;
      int code = 0;
      ret = pmodel->GetSpeciesSupportAndCode(sname, &supported, &code);
      if (ret)
        throw std::runtime_error("Error reading species support " + sname.ToString() + "(");
      if (supported)
        species.push_back(sname.ToString());
    }

    print_item("NUM_SPECIES", species.size());
    std::sort(species.begin(), species.end());
    print_item("SPECIES", species.begin(), species.end());
  }

  void print_cutoffs(kim_model_p& pmodel)
  {
    double influence_distance = 0.0;

    pmodel->GetInfluenceDistance(&influence_distance);

    print_item("INFLUENCE_DISTANCE", influence_distance);

    int num_lists = 0;
    const double *cutoffs_raw = nullptr;
    const int *data = nullptr;

    pmodel->GetNeighborListPointers(&num_lists, &cutoffs_raw, &data);

    print_item("NEIGHBOR_LIST_COUNT", num_lists);

    std::vector<double> cutoffs(cutoffs_raw, cutoffs_raw + num_lists);

    print_item("CUTOFFS", cutoffs.begin(), cutoffs.end());
  }

  void print_routines(kim_model_p& pmodel)
  {
    int present = 0;
    int required = 0;

    auto ret = pmodel->IsRoutinePresent(KIM::MODEL_ROUTINE_NAME::WriteParameterizedModel, &present, &required);
    if (ret)
      throw std::runtime_error("Error reading routines present");

    print_item("WRITE_MODEL", present);
  }

  void check_kim_model(const std::string& model_name)
  {
    kim_model_p pmodel(create_kim_model(model_name), [](KIM::Model* m){ KIM::Model::Destroy(&m); });

    int numberOfParameters = 0;

    pmodel->GetNumberOfParameters(&numberOfParameters);

    print_item("PARAMS", numberOfParameters);

    print_parameters(numberOfParameters, pmodel);

    print_species(pmodel);

    print_cutoffs(pmodel);

    print_routines(pmodel);
  }
}

int main(int argc, char* argv[]) {
  if (argc != 2 || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
    std::cout << g_usage << "\n";
    return EXIT_FAILURE;
  }

  try {
    check_kim_model(argv[1]);
  } catch (std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
