#ifndef __SYMMETRY_ANALYZER_DATABASE_HPP__
#define __SYMMETRY_ANALYZER_DATABASE_HPP__

#include <string>
#include <unordered_map>
#include <variant>
using hall_header_type = std::variant<bool, int, int, std::string>;
#include "csvreader/csv.h"

const std::vector<short> spglib_hall_numbers = {
    1,	 2,   3,   6,	9,   18,  21,  30,  39,	 57,  60,  63,	72,  81,  90,
    108, 109, 112, 115, 116, 119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
    155, 161, 164, 170, 173, 176, 182, 185, 191, 197, 203, 209, 212, 215, 218,
    221, 227, 228, 230, 233, 239, 245, 251, 257, 263, 266, 269, 275, 278, 284,
    290, 292, 298, 304, 310, 313, 316, 322, 334, 335, 337, 338, 341, 343, 349,
    350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 361, 363, 364, 366, 367,
    368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382,
    383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
    398, 399, 400, 401, 402, 404, 406, 407, 408, 410, 412, 413, 414, 416, 418,
    419, 420, 422, 424, 425, 426, 428, 430, 431, 432, 433, 435, 436, 438, 439,
    440, 441, 442, 443, 444, 446, 447, 448, 449, 450, 452, 454, 455, 456, 457,
    458, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
    475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489,
    490, 491, 492, 493, 494, 495, 497, 498, 500, 501, 502, 503, 504, 505, 506,
    507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 520, 521, 523,
    524, 525, 527, 529, 530};

class Hall_Table {
 public:
  enum class HALL_TABLE_HEADERS : int {
    HALL = 0,
    SPG_NUM = 1,
    SPG_FULL = 2,
    SYMBOL = 3,
    P = 4,
    Pinv = 5,
    PERM = 6
  };

  hall_header_type& operator[](HALL_TABLE_HEADERS header) {
    switch (header) {
      case HALL_TABLE_HEADERS::HALL:
	return std::ref(hall_numbers);
      case HALL_TABLE_HEADERS::SPG_NUM:
	return std::ref(spg_nums);
      case HALL_TABLE_HEADERS::SPG_FULL:
	return std::ref(spg_nums_full);
    }
  }

 private:
  std::vector<int> hall_numbers;
  std::vector<int> spg_nums;
  std::vector<std::string> spg_nums_full;
  std::vector<std::string> symbols;
  std::vector<std::string> Ps;
  std::vector<std::string> Pinvs;
  std::vector<bool> permutations;
};

CSVReader::io::CSVReader<7> hall_table("database/HM_full.csv");

#endif
