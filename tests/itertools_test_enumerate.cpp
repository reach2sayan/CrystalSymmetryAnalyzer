#include <catch2/catch_test_macros.hpp>

#include "base.hpp"
#include "enumerate.hpp"
#include "testtools.hpp"

namespace it = IteratorTools;
namespace tt = TestTools;

using Vec = std::vector<std::pair<std::size_t, char>>;
TEST_CASE("Basic Functioning enumerate", "[enumerate]") {
  std::vector<int> arr = {4, 5, 6};
  auto e = it::enumerate(arr);
  Vec v(std::begin(e), std::end(e));
  Vec vc{{0, 4}, {1, 5}, {2, 6}};

  REQUIRE(v == vc);
}
TEST_CASE("enumerate: can modify underlying sequence", "[enumerate]") {
  std::vector<char> s = {'a', 'b', 'c'};
  for (auto&& [i, c] : it::enumerate(s)) {
    c = '-';
  }
  REQUIRE(s[0] == '-');
  REQUIRE(s[1] == '-');
  REQUIRE(s[2] == '-');
}

TEST_CASE("enumerate: has .index, .element, .first, and .second") {
  std::vector<int> s = {1, 2, 3};
  auto e = it::enumerate(s);
  auto it = std::begin(e);
  REQUIRE(it->index == it->first);
  REQUIRE(&it->element == &it->second);
}

TEST_CASE("Empty enumerate", "[enumerate]") {
  std::vector<int> emp = {};
  auto e = it::enumerate(emp);
  REQUIRE(std::begin(e) == std::end(e));
}

TEST_CASE("Postfix ++ enumerate", "[enumerate]") {
  std::vector<int> s{1, 2, 3};
  auto e = it::enumerate(s);
  auto it = std::begin(e);
  it++;
  REQUIRE((*it).first == 1);
}

TEST_CASE("enumerate: structured bindings", "[enumerate]") {
  {
    std::string s{"amz"};
    auto e = it::enumerate(s);
    auto it = std::begin(e);
    REQUIRE(std::tuple_size<std::decay_t<decltype(*it)>>{} == 2);
    REQUIRE(std::get<0>(*it) == it->first);
  }

  Vec v;
  for (auto&& [i, c] : it::enumerate(std::string{"xyz"})) {
    v.emplace_back(i, c);
  }
  const Vec vc{{0, 'x'}, {1, 'y'}, {2, 'z'}};
  REQUIRE(v == vc);
}

TEST_CASE("Modifications through enumerate affect container", "[enumerate]") {
  std::vector<int> v{1, 2, 3, 4};
  std::vector<int> vc(v.size(), -1);
  for (auto&& p : it::enumerate(v)) {
    p.second = -1;
  }

  REQUIRE(v[0] == vc[0]);
}
TEST_CASE("enumerate with static array works", "[enumerate]") {
  std::vector<char> arr = {'w', 'x', 'y'};

  SECTION("Conversion to vector") {
    auto e = it::enumerate(arr);
    Vec v(std::begin(e), std::end(e));
    Vec vc{{0, 'w'}, {1, 'x'}, {2, 'y'}};
    REQUIRE(v == vc);
  }

  SECTION("Modification through enumerate") {
    for (auto&& p : it::enumerate(arr)) {
      p.second = 'z';
    }
    std::vector<char> v(std::begin(arr), std::end(arr));
    decltype(v) vc(v.size(), 'z');
    REQUIRE(v == vc);
  }
}

TEST_CASE("const enumerate", "[enumerate][const]") {
  Vec v;
  SECTION("lvalue") {
    std::string str = "abc";
    const auto e = it::enumerate(str);
    v.assign(std::begin(e), std::end(e));
  }
  SECTION("rvalue") {
    const auto e = it::enumerate(std::string("abc"));
    v.assign(std::begin(e), std::end(e));
  }
  SECTION("const lvalue") {
    const std::string str = "abc";
    const auto e = it::enumerate(str);
    v.assign(std::begin(e), std::end(e));
  }

  Vec vc{{0, 'a'}, {1, 'b'}, {2, 'c'}};

  REQUIRE(v == vc);
}
TEST_CASE("enumerate: operator->", "[enumerate]") {
  std::vector<int> ns = {50, 60, 70};
  auto e = it::enumerate(ns);
  auto it = std::begin(e);
  REQUIRE(it->first == 0);
  REQUIRE(it->second == 50);
}

TEST_CASE("enumerate: index and element", "[enumerate]") {
  std::string s{"ace"};
  auto e = it::enumerate(s);
  auto it = std::begin(e);
  REQUIRE((*it).index == 0);
  REQUIRE((*it).element == 'a');

  Vec v;
  for (auto&& p : it::enumerate(s)) {
    v.emplace_back(p.index, p.element);
  }
  Vec vc{{0, 'a'}, {1, 'c'}, {2, 'e'}};
  REQUIRE(v == vc);
}

TEST_CASE("enumerate: index and element through arrow", "[enumerate]") {
  std::string s{"ace"};
  auto e = it::enumerate(s);
  SECTION("One inspection") {
    auto it = std::begin(e);
    REQUIRE(it->index == 0);
    REQUIRE(it->element == 'a');
  }

  SECTION("full loop") {
    Vec v;
    for (auto it = std::begin(e), end_it = std::end(e); it != end_it; ++it) {
      v.emplace_back(it->index, it->element);
    }
    Vec vc{{0, 'a'}, {1, 'c'}, {2, 'e'}};
    REQUIRE(v == vc);
  }
}

TEST_CASE("Works with const iterable", "[enumerate]") {
  const std::string s{"ace"};
  auto e = it::enumerate(s);
  Vec v(std::begin(e), std::end(e));
  Vec vc{{0, 'a'}, {1, 'c'}, {2, 'e'}};
  REQUIRE(v == vc);
}

TEST_CASE("binds reference when it should", "[enumerate]") {
  tt::IterableType<char> bi{'x', 'y', 'z'};
  auto e = it::enumerate(bi);
  (void)e;
  REQUIRE_FALSE(bi.was_moved_from());
}

TEST_CASE("moves rvalues into enumerable object", "[enumerate]") {
  tt::IterableType<char> bi{'x', 'y', 'z'};
  // std::vector<char> bi{'x', 'y', 'z'};
  auto e = it::enumerate(std::move(bi));
  REQUIRE(bi.was_moved_from());
  (void)e;
}

TEST_CASE("Doesn't move or copy elements of iterable", "[enumerate]") {
  constexpr tt::MonolithObject<int> arr[] = {{6}, {7}, {8}};
  for (auto&& i : it::enumerate(arr)) {
    (void)i;
  }
}

TEST_CASE("enumerate: iterator meets requirements", "[enumerate]") {
  std::string s{};
  auto c = it::enumerate(s);
  REQUIRE(it::is_iterator_v<decltype(std::begin(c))>);
  REQUIRE(tt::reference_t_matches_deref_t<decltype(std::begin(c))>::value);
}

template <typename T>
using Imp_t = decltype(it::enumerate(std::declval<T>()));
TEST_CASE("enumerate: has correct ctor and assign ops", "[enumerate]") {
  REQUIRE(it::is_move_constructible_only<Imp_t<std::string&>>::value);
  REQUIRE(it::is_move_constructible_only<Imp_t<std::string>>::value);
}

