# Commercial launch drafts

Prepared 2026-10-16. **DRAFT: not approved, published, or binding.**
These proposed terms need owner review and appropriate legal/tax advice.
Do not copy unresolved fields into a live listing. No seller registration,
professional qualification, scientific certification, or sales is asserted.

## 1. Product and pricing proposal

Proposed public brand: Mongelli Labs, subject to trade-name verification.
The contracting seller's verified legal identity must be supplied through
the appropriate platform fields before publication.

| Item | Proposal | Release condition |
|---|---|---|
| Current alpha | Free technical evaluation; no paid launch now | Owner approves any public release; prominently retain limitations |
| Later educational bundle | USD 19 one-time starting-price hypothesis | Scientific scope, provenance, installation and commercial gates pass |
| Future OLED solver | Not priced or offered yet | Separate implementation, validation, scope and support review |
| Consulting/custom development | Separate written quote | No engineering certification, outcome or delivery date implied |

USD 19 is a positioning experiment, not market-researched willingness to pay.
Charge for a documented, convenient bundle or expressly agreed services,
not exclusive access to existing Apache-2.0 source. Do not imply a purchase
is required to exercise rights already granted by the open-source licenses.
Do not promise the future OLED solver as part of this purchase.

### Fee illustration, not net-income advice

Official Gumroad fees checked on 2026-10-16:

- Direct website sale: 10% + USD 0.50, with credit-card processing separately
  stated as 2.9% + USD 0.30. PayPal fees vary.
- Discover marketplace sale: flat 30%, including processing.
- Different account-specific rates, applicable taxes, currency/payout charges,
  refunds, and future fee changes can alter proceeds.

For a USD 19 direct card sale, the illustrative remainder is
`19 - (19 * 0.10 + 0.50) - (19 * 0.029 + 0.30) = USD 15.749`,
approximately USD 15.75 before other applicable adjustments and business taxes.
At 30% Discover fees, the comparable amount is USD 13.30.
Neither figure is a guaranteed payout or profit. Recheck live fees at launch.

Gumroad's merchant-of-record sales-tax handling does not eliminate the seller's
income/self-employment taxes, recordkeeping or other business obligations.

## 2. Storefront copy

**Title:** Mongelli Plane-Wave Optics Educational Preview

**Summary:** Explore thin-film transmission with a Python command-line tool.
Export numerical results, a contour plot, and reproducibility metadata.

**For:** Researchers, educators, and scientific programmers learning planar
multilayer optics. This is not a design-certified device simulator.

**Included in the proposed bundle:**

- Versioned Python wheel and source archive with license notices.
- Installation guide and command-line examples.
- Example PNG, CSV and JSON outputs.
- Validation notes covering specified plane-wave reference cases.
- Release checksums and runtime-dependency CycloneDX inventory.

**Requirements:** Python 3.11-3.14, ability to install Python packages, and
internet access to resolve dependencies unless separately supplied. Tests have
passed on GitHub-hosted Windows and Linux across that Python range.
macOS has not been tested. No standalone Windows installer or GUI is supplied.

**Model:** Air / illustrative organic film / 50 nm illustrative silver / air.
Default wavelengths: 400-700 nm; incidence angles: 0-89 degrees.
Transmittance is the equal-weight average of s and p polarizations.

**Limitations:** Alpha research/educational software. The index formulas are
illustrative and have no verified measurement provenance. No emitting dipoles,
calibrated OLED outcoupling, external quantum efficiency, device lifetime,
color accuracy, or fabrication-performance guarantee. Historical notebooks
are not part of the installable wheel. Do not describe the preview as the
future OLED solver.

**Price/support/refunds:** Insert only the owner-approved price and policies.
Keep this product unpublished until every required gate is complete.

## 3. Proposed customer terms

### Software rights and intended use

Project code remains Apache-2.0; project documentation and figures remain
CC BY 4.0 under the existing repository notices. Third-party components retain
their respective licenses. Preserve applicable attribution and notices.
The purchase, if offered, is nonexclusive and does not revoke or replace those
licenses. The educational description states the intended scope; it does not
add a field-of-use restriction to Apache-2.0.

No independently certified engineering result or device-performance warranty
is offered. Existing license disclaimers apply to the extent permitted by law;
mandatory consumer rights are not excluded. No additional bespoke liability
cap, arbitration clause, or governing-law clause is proposed here.

### Proposed support and updates policy

- Private support: use the seller contact shown on the Gumroad receipt.
  Confirm the mailbox and its monitoring before publishing.
- Public, non-sensitive reproducible bugs: use the repository issue tracker.
- Proposed response target: within three business days, Monday-Friday,
  excluding seller-observed US holidays. This is a target, not a guaranteed SLA.
- Proposed installation assistance window: 30 calendar days from purchase,
  limited to documented commands on supported platforms.
- No remote machine administration, confidential-data analysis, custom modeling,
  or scientific certification is included.
- Purchase delivers the stated version. Future fixes may be released, but
  perpetual support, all future versions, or the OLED solver are not promised.
- Announce critical known defects through available release/advisory channels;
  do not promise uninterrupted service or defect-free results.

These are proposed commitments only. The owner must confirm support capacity
before activation, including arrangements for absences or discontinuation.

### Proposed refund policy

Suggested customer-facing text, effective only after owner approval:

> Contact the seller through your Gumroad receipt within 30 calendar days of
> purchase to request a refund if the download cannot be used on a documented
> supported configuration or the product materially differs from its listing.
> Please provide your receipt reference, operating system, Python version and
> a non-sensitive error description. Troubleshooting is optional; we will not
> require confidential files. Refund requests are handled through Gumroad.
> Mandatory statutory rights and Gumroad's applicable policies remain in effect,
> including any legally required withdrawal rights.

Do not require deletion of Apache-licensed code or surrender of existing
license rights as a refund condition. Refund timing and payment-processing
charges are governed by the platform and payment rail, not guaranteed here.
Confirm the exact policy in the product's refund settings before launch.

### Proposed privacy notice

Draft operational commitments, not a statement that a live service exists:

> This product runs locally. The packaged preview has no application-level
> telemetry, cloud account, or license-server requirement. Installation may
> contact package repositories to retrieve dependencies. If you purchase through
> Gumroad, Gumroad and its payment providers process checkout, identity, payment
> and tax information under their own policies.
>
> The seller uses order details made available by Gumroad and information you
> voluntarily send for delivery, support, refunds and legally required records.
> Do not send tax IDs, banking information, identity documents, passwords or
> confidential research data to support. We do not need them to troubleshoot
> this local program. Do not post private data in public GitHub issues.
>
> Proposed support-ticket retention is 90 days after closure, except records
> needed for unresolved disputes, security incidents or legal obligations.
> Transaction records are retained as legally required. Contact the seller
> listed on your receipt for privacy requests; identity verification may be
> necessary. No marketing enrollment is required for product delivery.

Before use, verify actual support tools, providers, storage/access controls,
retention/deletion processes, public contact details, applicable rights and
jurisdictions. Do not claim deletion capabilities or compliance certifications
that have not been implemented. Gumroad's own retention is outside seller control.

## 4. Validation and provenance statement

Evidence for commit `58c8ab3dd7aff213e34886d018d828a6dfe0dce6`:

- [CI run](https://github.com/gmongell/Mathematica_python_TMM_OLEDPixel/actions/runs/38967043233):
  eight Windows/Linux Python 3.11-3.14 environments, 56 cases each,
  zero failures/errors/skips in parsed JUnit reports.
- Analytic Fresnel interface and Airy film checks, normal-incidence absorbing
  characteristic-matrix cross-check, lossless power balance and CLI tests.
- Wheel/source builds, metadata checks, isolated buyer installation, checksums
  and runtime dependency audit passed. An advisory scan is not a full security
  review and becomes stale as dependencies and advisories change.
- The optional draft-release job was not exercised. No public release is
  certified by this document.

| Evidence or asset | Current status | Remaining action |
|---|---|---|
| Original scripts and contributor rights | Existing repository notices; no complete ownership investigation | Owner reviews history, contributions and any institutional/employer obligations |
| Organic/silver formulas | Inherited placeholders, not material measurements | Trace origin or replace with clearly licensed measured data |
| `tmm`, NumPy, Matplotlib and transitive dependencies | Independently installed; inventory available | Review upstream license/NOTICE obligations for the precise delivery method |
| Scientific reference checks | Bounded plane-wave cases tested | Independent expert review and broader absorbing/oblique benchmarks |
| Images/screenshots | Locally generated preview output available | Select exact-version images with limitations visible |
| Documentation license | Existing CC BY 4.0 notice preserved | Confirm attribution and any adapted third-party content |
| Paid offering | No live listing | Seller verification, approved terms, delivery and checkout tests |

Existing tests do not establish measured-material accuracy or OLED emission
validity. Do not reuse these results as validation of newly added solver physics.

## 5. Owner approval record (all pending)

- [ ] Verify seller legal identity, trade name and authority to sell.
- [ ] Resolve code/contributor/material provenance and license review.
- [ ] Approve limited scientific claims after independent review.
- [ ] Approve price, delivery contents and screenshot set.
- [ ] Approve support capacity, contact, refund and updates policies.
- [ ] Review privacy/legal text against actual operations.
- [ ] Complete seller verification and accept agreements personally.
- [ ] Execute supported test checkout, verify receipt and buyer download.
- [ ] Approve exact release commit, tag/signing process and publication.

## Official references

Checked 2026-10-16; platform rules and account-specific requirements can change.

- [Fees](https://gumroad.com/help/article/66-gumroads-fees)
- [Pricing](https://gumroad.com/pricing)
- [Custom refund policy](https://gumroad.com/help/article/335-custom-refund-policy)
- [Testing a purchase](https://gumroad.com/help/article/62-testing-a-purchase)
- [Terms of service](https://gumroad.com/terms)
- [Privacy policy](https://gumroad.com/privacy)
