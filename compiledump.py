def func(x0, x1, tensor):
  temp_21584 = Interval(mpf('5.00000000000000000000e-1'), mpf('5.00000000000000000000e-1'))
  temp_21585 = 1.0
  temp_21586 = Interval(mpf('2.00000000000000000000e0'), mpf('2.00000000000000000000e0'))
  temp_21587 = Interval(mpf('-1.00000000000000000000e0'), mpf('-1.00000000000000000000e0'))
  temp_21588 = x0
  temp_21589 = (temp_21588 ** 2)
  temp_21590 = (temp_21587 * temp_21589)
  temp_21591 = x1
  temp_21592 = (temp_21591 ** 2)
  temp_21593 = (temp_21589 + temp_21592)
  temp_21594 = (temp_21593 ** 2)
  temp_21595 = (temp_21586 * temp_21594)
  temp_21596 = (temp_21590 / temp_21595)
  temp_21597 = temp_21593.sqrt()
  temp_21598 = (temp_21587 / temp_21597)
  temp_21599 = (temp_21596 - temp_21598)
  temp_21600 = (temp_21586 * temp_21599)
  temp_21601 = temp_21600.exp()
  temp_21602 = (temp_21585 / temp_21601)
  temp_21603 = (temp_21585 * temp_21601)
  temp_21604 = (temp_21585 * temp_21588)
  temp_21605 = (temp_21604 * temp_21585)
  temp_21606 = (temp_21586 * temp_21605)
  temp_21607 = (temp_21587 * temp_21606)
  temp_21608 = (temp_21585 * temp_21593)
  temp_21609 = (temp_21608 * temp_21606)
  temp_21610 = (temp_21586 * temp_21609)
  temp_21611 = (temp_21586 * temp_21610)
  temp_21612 = (temp_21611 * temp_21596)
  temp_21613 = (temp_21607 - temp_21612)
  temp_21614 = (temp_21613 / temp_21595)
  temp_21615 = (temp_21606 / temp_21586)
  temp_21616 = (temp_21615 / temp_21597)
  temp_21617 = (temp_21616 * temp_21598)
  temp_21618 = (-(temp_21617))
  temp_21619 = (temp_21618 / temp_21597)
  temp_21620 = (temp_21614 - temp_21619)
  temp_21621 = (temp_21586 * temp_21620)
  temp_21622 = (temp_21603 * temp_21621)
  temp_21623 = (temp_21622 + temp_21622)
  temp_21624 = (temp_21623 - temp_21622)
  temp_21625 = (temp_21602 * temp_21624)
  temp_21626 = (temp_21584 * temp_21625)
  tensor[(0, 0, 0)] = temp_21626
  temp_21627 = (temp_21585 * temp_21591)
  temp_21628 = (temp_21627 * temp_21585)
  temp_21629 = (temp_21586 * temp_21628)
  temp_21630 = (temp_21608 * temp_21629)
  temp_21631 = (temp_21586 * temp_21630)
  temp_21632 = (temp_21586 * temp_21631)
  temp_21633 = (temp_21632 * temp_21596)
  temp_21634 = (-(temp_21633))
  temp_21635 = (temp_21634 / temp_21595)
  temp_21636 = (temp_21629 / temp_21586)
  temp_21637 = (temp_21636 / temp_21597)
  temp_21638 = (temp_21637 * temp_21598)
  temp_21639 = (-(temp_21638))
  temp_21640 = (temp_21639 / temp_21597)
  temp_21641 = (temp_21635 - temp_21640)
  temp_21642 = (temp_21586 * temp_21641)
  temp_21643 = (temp_21603 * temp_21642)
  temp_21644 = (temp_21602 * temp_21643)
  temp_21645 = (temp_21584 * temp_21644)
  tensor[(0, 0, 1)] = temp_21645
  tensor[(0, 1, 0)] = temp_21645
  temp_21646 = (-(temp_21622))
  temp_21647 = (temp_21602 * temp_21646)
  temp_21648 = (temp_21584 * temp_21647)
  tensor[(0, 1, 1)] = temp_21648
  temp_21649 = (-(temp_21643))
  temp_21650 = (temp_21602 * temp_21649)
  temp_21651 = (temp_21584 * temp_21650)
  tensor[(1, 0, 0)] = temp_21651
  temp_21652 = (temp_21602 * temp_21622)
  temp_21653 = (temp_21584 * temp_21652)
  tensor[(1, 0, 1)] = temp_21653
  tensor[(1, 1, 0)] = temp_21653
  temp_21654 = (temp_21643 + temp_21643)
  temp_21655 = (temp_21654 - temp_21643)
  temp_21656 = (temp_21602 * temp_21655)
  temp_21657 = (temp_21584 * temp_21656)
  tensor[(1, 1, 1)] = temp_21657
  return
