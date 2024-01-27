In order to take into account the credit risk associated with lending or borrowing money in the interbank market at different time horizons, we focus our attention on a multicurve framework, which is
characterized by a discounting curve (i.e. a curve of discount factors that are used to actualize future cash flows) and several pseudo-discounting curves (one for each Euribor tenor), which lets us deal
with floating rates in financial contracts. The underlying idea of this choice is the large basis spreads
among different Euribor tenors, which makes problematic the use of a standard bootstrap (with only a
discounting curve) and introduces the need of a new bootstrap (section 3).
Since the aim of our analysis is the pricing of a contract in which fixed and floating rates are exchanged with quarterly payments in both legs, we are only interested in Euribor 3m. For this reason,
our setting presents a discounting curve and only one pseudo-discounting curve (related to Euribor
3m), leading us to a dual curve framework.
The contract that we aim to price is a so-called ”Amortizing Swap” between a bank and a corporate:
it is an interest rate derivative in which typically the payments of the fixed leg refer to a decreasing
notional, while the ones in the floating leg are all related to the initial value of the notional. However,
there exist also amortizing swaps where both the legs are characterized by a decreasing notional. We
denote the first case with ”single amortizing”, while the second with ”double amortizing”. Moreover,
these products are used in hedging instruments with declining principal, such as mortgages.
We first proceed with the risk free pricing of our contract (section 4) and then move on to a valuation
which embeds the CVA (section 5), which requires swaption pricing. In order to complete this task, we
need to utilize a suitable model for this multicurve framework. Since market data contain swaptions
quoted with Bachelier implied volatilities, we apply Bachelier closed formula, properly generalized to
this setting.
In section 6, we examine the case of an unwinding of the contract by the corporate: the aim is to
find the residual NPV of the instrument in that situation.
In order to employ an alternative method for swaption pricing, we have to consider that this particular context for discounting curve requires the structuring of an adequate model: we need to find
a parsimonious multicurve extension of Heath-Jarrow-Morton (HJM) models in order to adapt their
properties to this framework. In particular, we select the idea proposed in Baviera (2019), namely a
multicurve version of the Hull-White (MHW) model, characterized by only 3 parameters, which grants the possibility to avoid too complex calibrations and obtain simple closed formulas when pricing
interest rate derivatives (section 7).
Finally, we apply the closed formula coming from the MHW model calibrated in section 6 (obtained
through a generalization of the Jamshidian approach) in order to value the amortizing swap (with
CVA), as well as the trinomial tree, namely a numerical technique, based on the calibrated model itself
(section 8). Furthermore, in the same section we propose an alternative use of the Jamshidian trick for
the pricing formula.
The code for our analysis is realized in Matlab and Python.
