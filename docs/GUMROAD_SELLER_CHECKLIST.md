# Gumroad seller verification: guided checklist

Prepared 2026-10-16. This is operational guidance, not legal/tax advice or
confirmation of approval. Do not put completed identity, tax, address, banking
or document fields into this repository, chat, screenshots or GitHub issues.

## Current access and next action

Opening the official [payout settings](https://gumroad.com/settings/payments)
in the available browser redirected to Gumroad login on 2026-10-16.
No authenticated seller profile, verification status or product was accessed.
No form was submitted, account created, agreement accepted or purchase made.

The next user-only action is to sign in (or register if needed) directly with
Gumroad in a trusted browser, complete any MFA, and review the requested
onboarding steps. Do not send the assistant your password, one-time code,
recovery code, identity document, full account number, SSN or EIN.
Authentication in your usual browser may not authenticate the separate browser
available to this session.

### Preferred browser: existing Microsoft Edge session

At the owner's request, the official payout-settings URL was opened using the
running Microsoft Edge installation on 2026-10-16. The tool connection could
not observe/control those existing Edge tabs or confirm their login state.
No Edge process was terminated, profile copied, debugging mode enabled,
extension accessed, or saved password extracted.

Use Proton Pass normally inside Edge if sign-in is requested. Complete MFA,
identity/tax/bank entry and agreement acceptance yourself on the authenticated
official page. The assistant can interpret non-sensitive instructions/status,
but cannot infer approval or submit invisible forms.

For help afterward, provide only the non-sensitive field label or error text
and a status such as "identity pending", "bank requested", or "account review".
Omit numbers, personal document images and private dashboard URLs/tokens.

## 1. Select the truthful account type

Gumroad's current guidance says:

- **Individual:** operating under your name, effectively as a sole proprietor,
  or without business registration documents.
- **Business:** you have business registration documents; a representative's
  personal details and registration documents may also be required.

The earlier enrollment packets described a sole proprietorship, but this is not
independent evidence of its present legal or registration status. Do not select
"company", claim incorporation, or use a trade name as a substitute for the
verified legal payee unless the actual portal and documents support it.
Use Gumroad support for account-type ambiguity rather than inventing an entity.

## 2. Prepare privately and enter only in the official portal

| Portal requirement | Preparation | Who completes it |
|---|---|---|
| Country | Actual supported payout country | Account owner |
| Legal name | Exact name, including middle name where required, matching official records | Account owner |
| Physical address | Valid actual address; Gumroad excludes PO boxes and military addresses in its guidance | Account owner |
| Date of birth/phone | Accurate details; complete any requested verification | Account owner |
| Tax identity | Name/TIN matching IRS records; individual/sole-proprietor treatment depends on actual facts | Account owner, with tax adviser if needed |
| Business evidence | Registration and representative/ownership evidence only if applicable/requested | Authorized representative |
| Photo ID/address evidence | Valid, readable documentation if requested by the secure verification flow | Account owner |
| Payout account | Actual eligible account and payee details matching the selected type | Account owner |
| Agreements | Read and accept personally only with authority | Account owner |

For US sellers, Gumroad describes EIN use and SSN use for individuals/sole
proprietors without an EIN. Follow the live field instructions and IRS W-9
instructions for the correct name/TIN combination; do not infer it from a brand.
An EIN does not by itself prove incorporation or remove all personal KYC.

Gumroad uses Stripe for card processing/bank payouts; the secure flow can request
additional identity verification. PayPal availability depends on the account
and country, not preference alone. Do not switch countries to obtain a payout
method. Gumroad warns country changes can clear saved details and forfeit some
unpaid amounts; contact support first if your country is wrong.

## 3. Verify each separate status

- [ ] Login/email and account security complete.
- [ ] Account type and country accurately selected.
- [ ] Legal identity submitted and verification status checked.
- [ ] Tax information submitted and accepted as required.
- [ ] Payout account added and verification status checked.
- [ ] Any business/representative evidence accepted.
- [ ] Any outstanding dashboard warnings resolved.
- [ ] Seller/account review status understood.
- [ ] Product publication eligibility confirmed separately from payouts.

Do not mark these checks merely because a form was saved. Identity verification,
product eligibility, fraud/account review, and funds actually arriving are
different milestones. Only the platform can confirm acceptance.

Gumroad's public payout guidance currently describes a standard USD 100 payout
minimum (higher in some locations), holding periods depending on payout schedule,
and possible account review taking 1-3 weeks. These are not guaranteed dates,
account-specific approval criteria, or a reason to create artificial transactions.
Check your actual dashboard and current documentation.

## 4. Product setup after approval gates

Keep the product unpublished during preparation. Use
[commercial drafts](COMMERCIAL_DRAFTS.md) as proposed copy, not already approved
terms. Confirm the real product type, supported platforms, price, support
contact, license notices, refund settings and privacy link.

Before uploading:

1. Obtain the exact tested wheel/source from the approved build.
2. Verify checksums and audit inventory; include installation/license information.
3. Ensure no customer data, credentials, seller records or development environment
   is in the upload. Do not upload `.venv` or a whole workspace backup.
4. Show the alpha and material-data limitations before checkout.
5. Do not promise the not-yet-created OLED solver as a purchase entitlement.
6. Do not enable publication, paid subscriptions or recurring commitments merely
   to test a draft.

## 5. Supported checkout test

Official Gumroad instructions: remain logged in as the seller, open the product
purchase page and use the built-in **Test card** flow. If the test mode is not
visible, stop and use Gumroad's current supported instructions/support.
The help article also describes a 100% discount-code alternative.

- Never pay for your own product with your real card to trigger account review
  or test payouts. Gumroad explicitly warns this can cause suspension.
- If testing requires publishing or changing product visibility, obtain a
  deliberate publication decision first; do not silently make a draft public.
- Check price, tax display, refund-policy visibility, receipt/support link,
  filenames, licenses, checksums and download delivery.
- Install the buyer download in a clean environment and run the documented
  command. Confirm PNG/CSV/JSON contents against the shipped example.
- A test transaction verifies delivery, not actual settlement or KYC approval.

## 6. Safe support request template

> Subject: Seller onboarding/account-type clarification
>
> I am preparing a downloadable scientific/educational Python software product.
> I need help confirming the appropriate account type and the next onboarding
> step for my actual legal status and registration documents.
> The non-sensitive dashboard status is: [insert status only].
> Please direct me to the secure in-account process for any required evidence.
> I will not send tax identifiers, banking details or identity documents in
> ordinary email.

Send through Gumroad's official support channel from your authenticated account.
Do not claim production validation or a legal status that has not been verified.

## Official references

Checked 2026-10-16; the live portal controls actual requirements.

- [Payout settings guidance](https://gumroad.com/help/article/260-your-payout-settings-page)
- [Getting paid/account review](https://gumroad.com/help/article/13-getting-paid)
- [Test purchase](https://gumroad.com/help/article/62-testing-a-purchase)
- [Fees](https://gumroad.com/help/article/66-gumroads-fees)
- [Custom refund settings](https://gumroad.com/help/article/335-custom-refund-policy)
- [IRS W-9](https://www.irs.gov/forms-pubs/about-form-w-9)
